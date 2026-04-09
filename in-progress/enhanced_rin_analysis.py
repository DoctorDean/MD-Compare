#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import warnings

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, distance_array
    from MDAnalysis.analysis import contacts
except ImportError:
    print("Error: MDAnalysis is required.")
    sys.exit(1)

try:
    import networkx as nx
except ImportError:
    print("Error: networkx is required.")
    sys.exit(1)

try:
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import pearsonr
except ImportError:
    print("Warning: scipy not available. Some analyses will be limited.")

# =====================================================
# ENHANCED CONTACT MAP CALCULATION
# =====================================================


class EnhancedRIN:
    def __init__(self, topology, trajectory, selection="protein and not name H*"):
        """
        Enhanced Residue Interaction Network Analysis
        
        Parameters:
        -----------
        topology : str
            Path to topology file
        trajectory : str  
            Path to trajectory file
        selection : str
            MDAnalysis selection string for atoms to analyze
        """
        try:
            self.u = mda.Universe(topology, trajectory)
        except Exception as e:
            raise RuntimeError(f"Failed to load trajectory: {e}")
            
        self.atoms = self.u.select_atoms(selection)
        if len(self.atoms) == 0:
            raise RuntimeError("No atoms selected. Check selection string.")
            
        # Set up residue mapping
        self.resids = self.atoms.resids
        self.unique_res = np.unique(self.resids)
        self.n_res = len(self.unique_res)
        self.resid_to_index = {res: i for i, res in enumerate(self.unique_res)}
        self.atom_res_index = np.array([self.resid_to_index[r] for r in self.resids])
        
        print(f"Loaded system with {self.n_res} residues, {len(self.atoms)} atoms")
        print(f"Trajectory: {len(self.u.trajectory)} frames")

    def compute_multi_contact_maps(self, cutoffs={'all_atom': 4.5, 'ca_only': 8.0}, 
                                  interaction_types=['distance', 'hbond', 'salt_bridge']):
        """
        Compute multiple types of contact maps
        
        Parameters:
        -----------
        cutoffs : dict
            Different cutoff distances for different contact types
        interaction_types : list
            Types of interactions to analyze
        """
        self.contact_maps = {}
        self.interaction_details = {}
        
        # Initialize contact matrices
        for itype in interaction_types:
            self.contact_maps[itype] = np.zeros((self.n_res, self.n_res))
        
        n_frames = 0
        
        print("Computing contact maps...")
        for i, ts in enumerate(self.u.trajectory):
            if i % 100 == 0:
                print(f"  Frame {i}/{len(self.u.trajectory)}")
                
            # Distance-based contacts
            if 'distance' in interaction_types:
                self._compute_distance_contacts(cutoffs['all_atom'])
            
            # Hydrogen bonds
            if 'hbond' in interaction_types:
                self._compute_hydrogen_bonds()
                
            # Salt bridges  
            if 'salt_bridge' in interaction_types:
                self._compute_salt_bridges()
                
            n_frames += 1
        
        # Normalize by number of frames
        for itype in interaction_types:
            self.contact_maps[itype] /= n_frames
            
        print(f"Contact maps computed over {n_frames} frames")
        return self.contact_maps

    def _compute_distance_contacts(self, cutoff):
        """Compute distance-based contacts"""
        coords = self.atoms.positions
        
        pairs = capped_distance(
            coords, coords,
            max_cutoff=cutoff,
            return_distances=False
        )
        
        if len(pairs) > 0:
            atom_i, atom_j = pairs[:, 0], pairs[:, 1]
            res_i = self.atom_res_index[atom_i]
            res_j = self.atom_res_index[atom_j]
            
            # Only count inter-residue contacts
            mask = res_i != res_j
            res_i, res_j = res_i[mask], res_j[mask]
            
            self.contact_maps['distance'][res_i, res_j] += 1
            self.contact_maps['distance'][res_j, res_i] += 1

    def _compute_hydrogen_bonds(self, distance_cutoff=3.5, angle_cutoff=150):
        """Compute hydrogen bonds between residues"""
        # Simplified H-bond detection - can be enhanced with more sophisticated methods
        donors = self.u.select_atoms("protein and (name N or name O)")
        acceptors = self.u.select_atoms("protein and (name O or name N)")
        
        if len(donors) == 0 or len(acceptors) == 0:
            return
            
        # Distance-based H-bond screening
        donor_coords = donors.positions
        acceptor_coords = acceptors.positions
        
        distances = distance_array(donor_coords, acceptor_coords)
        hbond_pairs = np.where(distances <= distance_cutoff)
        
        for i, j in zip(hbond_pairs[0], hbond_pairs[1]):
            donor_res = self.resid_to_index.get(donors[i].resid)
            acceptor_res = self.resid_to_index.get(acceptors[j].resid)
            
            if donor_res is not None and acceptor_res is not None and donor_res != acceptor_res:
                self.contact_maps['hbond'][donor_res, acceptor_res] += 1
                self.contact_maps['hbond'][acceptor_res, donor_res] += 1

    def _compute_salt_bridges(self, cutoff=4.0):
        """Compute salt bridges between charged residues"""
        positive = self.u.select_atoms("protein and (resname ARG LYS) and (name NH* or name NZ)")
        negative = self.u.select_atoms("protein and (resname ASP GLU) and (name OD* or name OE*)")
        
        if len(positive) == 0 or len(negative) == 0:
            return
            
        pos_coords = positive.positions
        neg_coords = negative.positions
        
        distances = distance_array(pos_coords, neg_coords)
        salt_pairs = np.where(distances <= cutoff)
        
        for i, j in zip(salt_pairs[0], salt_pairs[1]):
            pos_res = self.resid_to_index.get(positive[i].resid)
            neg_res = self.resid_to_index.get(negative[j].resid)
            
            if pos_res is not None and neg_res is not None:
                self.contact_maps['salt_bridge'][pos_res, neg_res] += 1
                self.contact_maps['salt_bridge'][neg_res, pos_res] += 1

    def compute_advanced_network_metrics(self, interaction_type='distance', threshold=0.2):
        """
        Compute comprehensive network analysis metrics
        """
        contact_matrix = self.contact_maps[interaction_type]
        
        # Build network
        G = self._build_network(contact_matrix, threshold)
        
        metrics = {}
        
        # Basic network properties
        metrics['n_nodes'] = G.number_of_nodes()
        metrics['n_edges'] = G.number_of_edges()
        metrics['density'] = nx.density(G)
        
        # Centrality measures
        print("Computing centrality measures...")
        metrics['betweenness_centrality'] = nx.betweenness_centrality(G, weight='weight')
        metrics['closeness_centrality'] = nx.closeness_centrality(G, distance='weight')
        metrics['eigenvector_centrality'] = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
        metrics['degree_centrality'] = nx.degree_centrality(G)
        
        # Community detection
        print("Detecting communities...")
        try:
            communities = nx.community.greedy_modularity_communities(G, weight='weight')
            metrics['communities'] = [list(c) for c in communities]
            metrics['modularity'] = nx.community.modularity(G, communities, weight='weight')
        except:
            metrics['communities'] = []
            metrics['modularity'] = 0
        
        # Path analysis
        print("Computing path metrics...")
        if nx.is_connected(G):
            metrics['average_path_length'] = nx.average_shortest_path_length(G, weight='weight')
            metrics['diameter'] = nx.diameter(G)
        else:
            # For disconnected graphs
            largest_cc = max(nx.connected_components(G), key=len)
            G_cc = G.subgraph(largest_cc)
            metrics['average_path_length'] = nx.average_shortest_path_length(G_cc, weight='weight')
            metrics['diameter'] = nx.diameter(G_cc)
            
        # Clustering
        metrics['clustering_coefficient'] = nx.average_clustering(G, weight='weight')
        
        # Store network
        metrics['network'] = G
        
        return metrics

    def _build_network(self, contact_matrix, threshold):
        """Build NetworkX graph from contact matrix"""
        G = nx.Graph()
        
        # Add nodes
        for res in self.unique_res:
            G.add_node(int(res))
        
        # Add edges
        for i in range(self.n_res):
            for j in range(i + 1, self.n_res):
                weight = contact_matrix[i, j]
                if weight >= threshold:
                    G.add_edge(
                        int(self.unique_res[i]),
                        int(self.unique_res[j]),
                        weight=float(weight)
                    )
        
        return G

    def find_allosteric_pathways(self, source_residues, target_residues, 
                               interaction_type='distance', threshold=0.2, max_paths=5):
        """
        Find allosteric communication pathways between residue sets
        
        Parameters:
        -----------
        source_residues : list
            Source residues (e.g., mutation sites)
        target_residues : list
            Target residues (e.g., active site)
        interaction_type : str
            Type of interaction network to use
        threshold : float
            Contact frequency threshold for edges
        max_paths : int
            Maximum number of paths to find per source-target pair
        """
        contact_matrix = self.contact_maps[interaction_type]
        G = self._build_network(contact_matrix, threshold)
        
        pathways = {}
        
        for source in source_residues:
            for target in target_residues:
                if source in G and target in G:
                    try:
                        # Find shortest paths
                        paths = list(nx.shortest_simple_paths(G, source, target, weight='weight'))
                        pathways[f'{source}_to_{target}'] = {
                            'paths': paths[:max_paths],
                            'shortest_distance': nx.shortest_path_length(G, source, target, weight='weight')
                        }
                    except nx.NetworkXNoPath:
                        pathways[f'{source}_to_{target}'] = {
                            'paths': [],
                            'shortest_distance': float('inf')
                        }
                else:
                    pathways[f'{source}_to_{target}'] = {
                        'paths': [],
                        'shortest_distance': float('inf')
                    }
        
        return pathways

    def analyze_sub_trajectories(self, n_segments=5):
        """
        Analyze trajectory in segments for statistical robustness
        """
        total_frames = len(self.u.trajectory)
        segment_size = total_frames // n_segments
        
        segment_results = {}
        
        for seg in range(n_segments):
            print(f"Analyzing segment {seg+1}/{n_segments}")
            
            start_frame = seg * segment_size
            end_frame = start_frame + segment_size
            
            # Temporarily modify trajectory slice
            original_traj = self.u.trajectory
            segment_traj = self.u.trajectory[start_frame:end_frame]
            
            # Compute contact maps for this segment
            segment_contacts = {}
            for itype in ['distance']:  # Can extend to other types
                segment_contacts[itype] = np.zeros((self.n_res, self.n_res))
            
            n_frames = 0
            for ts in segment_traj:
                self._compute_distance_contacts(4.5)
                # Add to segment-specific matrix
                segment_contacts['distance'] += self.contact_maps['distance']
                n_frames += 1
            
            # Normalize
            for itype in segment_contacts:
                segment_contacts[itype] /= n_frames
            
            # Compute metrics for this segment
            temp_contact_maps = self.contact_maps
            self.contact_maps = segment_contacts
            
            segment_metrics = self.compute_advanced_network_metrics()
            segment_results[f'segment_{seg}'] = {
                'contact_maps': segment_contacts.copy(),
                'metrics': segment_metrics
            }
            
            # Restore original contact maps
            self.contact_maps = temp_contact_maps
        
        return segment_results

    def compute_residue_flexibility_correlation(self):
        """
        Compute correlation between network centrality and structural flexibility
        """
        # Get RMSF for CA atoms
        ca_atoms = self.u.select_atoms("name CA")
        rmsf_analysis = mda.analysis.rms.RMSF(ca_atoms)
        rmsf_analysis.run()
        
        rmsf_values = rmsf_analysis.rmsf
        
        # Get centrality values
        metrics = self.compute_advanced_network_metrics()
        centrality = metrics['betweenness_centrality']
        
        # Match residues
        correlations = {}
        for i, resid in enumerate(self.unique_res):
            if resid in ca_atoms.resids:
                ca_index = np.where(ca_atoms.resids == resid)[0][0]
                rmsf_val = rmsf_values[ca_index]
                cent_val = centrality.get(resid, 0)
                
                correlations[resid] = {
                    'rmsf': rmsf_val,
                    'centrality': cent_val
                }
        
        # Calculate overall correlation
        rmsf_vals = [v['rmsf'] for v in correlations.values()]
        cent_vals = [v['centrality'] for v in correlations.values()]
        
        if len(rmsf_vals) > 3:
            correlation_coeff, p_value = pearsonr(rmsf_vals, cent_vals)
        else:
            correlation_coeff, p_value = 0, 1
            
        return correlations, correlation_coeff, p_value

    def save_comprehensive_outputs(self, prefix, metrics=None):
        """
        Save all analysis outputs
        """
        # Save contact maps
        for itype, contact_matrix in self.contact_maps.items():
            np.save(f"{prefix}_{itype}_contact_map.npy", contact_matrix)
            np.savetxt(f"{prefix}_{itype}_contact_map.csv", contact_matrix, delimiter=",")
            
            # Enhanced visualization
            plt.figure(figsize=(10, 8))
            
            # Use better colormap and add residue labels
            im = plt.imshow(contact_matrix, origin="lower", cmap='viridis', 
                          extent=[self.unique_res.min(), self.unique_res.max(),
                                 self.unique_res.min(), self.unique_res.max()])
            
            plt.colorbar(im, label="Contact Frequency", shrink=0.8)
            plt.xlabel("Residue Number")
            plt.ylabel("Residue Number")
            plt.title(f"Residue Contact Map - {itype.replace('_', ' ').title()}")
            
            # Add grid for better readability
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"{prefix}_{itype}_contact_map.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        # Save network metrics if provided
        if metrics:
            # Save network
            nx.write_graphml(metrics['network'], f"{prefix}_network.graphml")
            
            # Save centrality measures
            centrality_data = []
            for res in self.unique_res:
                centrality_data.append([
                    res,
                    metrics['betweenness_centrality'].get(res, 0),
                    metrics['closeness_centrality'].get(res, 0),
                    metrics['eigenvector_centrality'].get(res, 0),
                    metrics['degree_centrality'].get(res, 0)
                ])
            
            np.savetxt(f"{prefix}_centrality_metrics.csv", 
                      centrality_data,
                      delimiter=",",
                      header="residue,betweenness,closeness,eigenvector,degree",
                      comments="")
            
            # Save network properties
            with open(f"{prefix}_network_properties.txt", 'w') as f:
                f.write(f"Number of nodes: {metrics['n_nodes']}\n")
                f.write(f"Number of edges: {metrics['n_edges']}\n")
                f.write(f"Network density: {metrics['density']:.4f}\n")
                f.write(f"Average clustering coefficient: {metrics['clustering_coefficient']:.4f}\n")
                f.write(f"Average path length: {metrics['average_path_length']:.4f}\n")
                f.write(f"Network diameter: {metrics['diameter']:.4f}\n")
                f.write(f"Modularity: {metrics['modularity']:.4f}\n")
                f.write(f"Number of communities: {len(metrics['communities'])}\n")

    def create_pymol_visualization_script(self, prefix, metrics):
        """
        Create PyMOL script for network visualization
        """
        script_content = f"""
# PyMOL script for network visualization
# Generated automatically

# Load structure
load {prefix}_topology.pdb

# Color by betweenness centrality
"""
        
        # Add coloring commands based on centrality
        centrality = metrics['betweenness_centrality']
        max_cent = max(centrality.values()) if centrality else 1
        
        for res, cent_val in centrality.items():
            normalized_cent = cent_val / max_cent
            color_val = int(normalized_cent * 9)  # 0-9 color scale
            script_content += f"color tv_blue, resi {res}\n"
            script_content += f"set_color cent_{res}, [{normalized_cent}, 0, {1-normalized_cent}]\n"
            script_content += f"color cent_{res}, resi {res}\n"
        
        script_content += """
# Show as cartoon
show cartoon
hide lines
show sticks, name CA
"""
        
        with open(f"{prefix}_pymol_network.pml", 'w') as f:
            f.write(script_content)

# =====================================================
# MAIN FUNCTION
# =====================================================

def main():
    parser = argparse.ArgumentParser(
        description="Enhanced residue interaction network analysis for MD trajectories."
    )
    
    parser.add_argument("-t", "--topology", required=True,
                       help="Topology file (PDB, GRO, etc.)")
    parser.add_argument("-x", "--trajectory", required=True,
                       help="Trajectory file (XTC, DCD, etc.)")
    parser.add_argument("-o", "--output", default="enhanced_rin",
                       help="Output prefix (default: enhanced_rin)")
    parser.add_argument("--cutoff", type=float, default=4.5,
                       help="Contact cutoff in Å (default: 4.5)")
    parser.add_argument("--threshold", type=float, default=0.2,
                       help="Contact frequency threshold for network edges (default: 0.2)")
    parser.add_argument("--selection", default="protein and not name H*",
                       help="Atom selection (default: protein and not name H*)")
    parser.add_argument("--segments", type=int, default=5,
                       help="Number of trajectory segments for statistical analysis (default: 5)")
    parser.add_argument("--source-residues", nargs='+', type=int,
                       help="Source residues for pathway analysis (e.g., mutation sites)")
    parser.add_argument("--target-residues", nargs='+', type=int,
                       help="Target residues for pathway analysis (e.g., active site)")
    parser.add_argument("--interaction-types", nargs='+', 
                       choices=['distance', 'hbond', 'salt_bridge'],
                       default=['distance'],
                       help="Types of interactions to analyze")
    
    args = parser.parse_args()
    
    try:
        # Initialize enhanced RIN analysis
        rin = EnhancedRIN(args.topology, args.trajectory, args.selection)
        
        # Compute contact maps
        cutoffs = {'all_atom': args.cutoff, 'ca_only': args.cutoff * 1.8}
        contact_maps = rin.compute_multi_contact_maps(cutoffs, args.interaction_types)
        
        # Compute network metrics
        print("Computing network metrics...")
        metrics = rin.compute_advanced_network_metrics('distance', args.threshold)
        
        # Trajectory segmentation analysis
        if args.segments > 1:
            print(f"Analyzing {args.segments} trajectory segments...")
            segment_results = rin.analyze_sub_trajectories(args.segments)
            
            # Save segment analysis
            import pickle
            with open(f"{args.output}_segment_analysis.pkl", 'wb') as f:
                pickle.dump(segment_results, f)
        
        # Pathway analysis
        if args.source_residues and args.target_residues:
            print("Finding allosteric pathways...")
            pathways = rin.find_allosteric_pathways(
                args.source_residues, args.target_residues
            )
            
            # Save pathway results
            with open(f"{args.output}_pathways.txt", 'w') as f:
                for pathway_name, pathway_data in pathways.items():
                    f.write(f"\n{pathway_name}:\n")
                    f.write(f"  Shortest distance: {pathway_data['shortest_distance']:.4f}\n")
                    for i, path in enumerate(pathway_data['paths'][:3]):
                        f.write(f"  Path {i+1}: {' -> '.join(map(str, path))}\n")
        
        # Flexibility correlation analysis
        print("Computing flexibility-centrality correlation...")
        try:
            flex_corr, corr_coeff, p_val = rin.compute_residue_flexibility_correlation()
            print(f"Flexibility-centrality correlation: r={corr_coeff:.3f}, p={p_val:.4f}")
        except Exception as e:
            print(f"Flexibility correlation failed: {e}")
        
        # Save all outputs
        print("Saving outputs...")
        rin.save_comprehensive_outputs(args.output, metrics)
        rin.create_pymol_visualization_script(args.output, metrics)
        
        print(f"\nAnalysis complete. Outputs saved with prefix: {args.output}")
        print(f"Network: {metrics['n_nodes']} nodes, {metrics['n_edges']} edges")
        print(f"Network density: {metrics['density']:.4f}")
        print(f"Average clustering: {metrics['clustering_coefficient']:.4f}")
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
