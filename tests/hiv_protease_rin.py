#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from enhanced_rin_analysis import EnhancedRIN
import argparse
import sys

class HIVProteaseRIN(EnhancedRIN):
    """
    HIV Protease-specific Residue Interaction Network Analysis
    """
    
    def __init__(self, topology, trajectory):
        super().__init__(topology, trajectory, selection="protein and not name H*")
        
        # HIV-1 protease specific residue definitions
        self.active_site = [25, 25]  # Asp25, Asp25' (assuming chain A and B)
        self.flap_regions = {
            'flap_tips': [50, 150],  # Ile50, Ile50' (chain A: 50, chain B: 150)  
            'flap_hinge': [37, 137, 38, 138, 39, 139, 40, 140]  # Hinge region
        }
        
        # Common resistance mutation sites
        self.primary_resistance_sites = [30, 32, 46, 47, 48, 50, 82, 84, 88, 90]
        self.secondary_resistance_sites = [10, 20, 24, 36, 53, 54, 60, 62, 71, 73, 77, 93]
        
        # Structural domains
        self.domains = {
            'N_terminal': list(range(1, 10)),
            'beta_hairpin_1': list(range(10, 23)),  
            'active_site_loop': list(range(23, 30)),
            'beta_hairpin_2': list(range(30, 43)),
            'flap': list(range(43, 58)),
            'alpha_helix': list(range(58, 75)),
            'beta_sheet': list(range(75, 85)),
            'C_terminal': list(range(85, 99))
        }
        
        print("HIV Protease RIN Analysis initialized")
        print(f"Active site residues: {self.active_site}")
        print(f"Flap regions defined: {list(self.flap_regions.keys())}")

    def compute_hiv_specific_metrics(self, interaction_type='distance', threshold=0.2):
        """
        Compute HIV protease-specific network metrics
        """
        # Get basic network metrics
        metrics = self.compute_advanced_network_metrics(interaction_type, threshold)
        
        # Add HIV-specific analyses
        hiv_metrics = {}
        
        # Flap-to-active site communication
        hiv_metrics['flap_active_pathways'] = self.find_allosteric_pathways(
            self.flap_regions['flap_tips'], 
            self.active_site,
            interaction_type, threshold
        )
        
        # Resistance site connectivity
        hiv_metrics['primary_resistance_centrality'] = {}
        hiv_metrics['secondary_resistance_centrality'] = {}
        
        for res in self.primary_resistance_sites:
            if res in metrics['betweenness_centrality']:
                hiv_metrics['primary_resistance_centrality'][res] = metrics['betweenness_centrality'][res]
        
        for res in self.secondary_resistance_sites:
            if res in metrics['betweenness_centrality']:
                hiv_metrics['secondary_resistance_centrality'][res] = metrics['betweenness_centrality'][res]
        
        # Inter-domain communication
        hiv_metrics['interdomain_communication'] = self._analyze_interdomain_communication(metrics['network'])
        
        # Flap dynamics network properties
        hiv_metrics['flap_network_properties'] = self._analyze_flap_network_properties(metrics['network'])
        
        return {**metrics, **hiv_metrics}

    def _analyze_interdomain_communication(self, network):
        """
        Analyze communication between structural domains
        """
        interdomain_edges = {}
        domain_centralities = {}
        
        # Calculate average centrality for each domain
        for domain_name, residues in self.domains.items():
            domain_nodes = [r for r in residues if r in network]
            if domain_nodes:
                centralities = [network.nodes[node].get('centrality', 0) for node in domain_nodes if node in network]
                domain_centralities[domain_name] = np.mean(centralities) if centralities else 0
        
        # Count inter-domain edges
        for domain1, residues1 in self.domains.items():
            for domain2, residues2 in self.domains.items():
                if domain1 < domain2:  # Avoid double counting
                    edge_count = 0
                    total_weight = 0
                    
                    for r1 in residues1:
                        for r2 in residues2:
                            if network.has_edge(r1, r2):
                                edge_count += 1
                                total_weight += network[r1][r2]['weight']
                    
                    interdomain_edges[f"{domain1}_{domain2}"] = {
                        'edge_count': edge_count,
                        'total_weight': total_weight,
                        'avg_weight': total_weight / edge_count if edge_count > 0 else 0
                    }
        
        return {
            'interdomain_edges': interdomain_edges,
            'domain_centralities': domain_centralities
        }

    def _analyze_flap_network_properties(self, network):
        """
        Analyze network properties specific to flap regions
        """
        flap_properties = {}
        
        # Combine all flap residues
        all_flap_residues = []
        for flap_res_list in self.flap_regions.values():
            all_flap_residues.extend(flap_res_list)
        
        flap_nodes = [r for r in all_flap_residues if r in network]
        
        if len(flap_nodes) > 0:
            # Flap subgraph
            flap_subgraph = network.subgraph(flap_nodes)
            
            flap_properties = {
                'flap_density': flap_subgraph.number_of_edges() / (len(flap_nodes) * (len(flap_nodes) - 1) / 2) if len(flap_nodes) > 1 else 0,
                'flap_clustering': np.mean([network.degree[node] for node in flap_nodes]) if flap_nodes else 0,
                'flap_to_nonflap_edges': sum(1 for node in flap_nodes for neighbor in network.neighbors(node) if neighbor not in flap_nodes)
            }
        
        return flap_properties

    def analyze_mutation_effects(self, mutation_sites, interaction_type='distance', threshold=0.2):
        """
        Analyze the network effects of specific mutations
        
        Parameters:
        -----------
        mutation_sites : list
            List of residue numbers where mutations occur
        """
        print(f"Analyzing mutation effects for sites: {mutation_sites}")
        
        # Get baseline network metrics
        baseline_metrics = self.compute_hiv_specific_metrics(interaction_type, threshold)
        baseline_network = baseline_metrics['network']
        
        mutation_analysis = {}
        
        for mut_site in mutation_sites:
            if mut_site not in baseline_network:
                print(f"Warning: Mutation site {mut_site} not in network")
                continue
                
            site_analysis = {}
            
            # Local network properties
            neighbors = list(baseline_network.neighbors(mut_site))
            site_analysis['local_neighbors'] = neighbors
            site_analysis['local_degree'] = baseline_network.degree[mut_site]
            site_analysis['local_clustering'] = nx.clustering(baseline_network, mut_site)
            
            # Centrality measures
            site_analysis['betweenness_centrality'] = baseline_metrics['betweenness_centrality'].get(mut_site, 0)
            site_analysis['closeness_centrality'] = baseline_metrics['closeness_centrality'].get(mut_site, 0)
            site_analysis['eigenvector_centrality'] = baseline_metrics['eigenvector_centrality'].get(mut_site, 0)
            
            # Distance to important sites
            site_analysis['distances_to_active_site'] = {}
            for active_res in self.active_site:
                if active_res in baseline_network:
                    try:
                        dist = nx.shortest_path_length(baseline_network, mut_site, active_res, weight='weight')
                        site_analysis['distances_to_active_site'][active_res] = dist
                    except nx.NetworkXNoPath:
                        site_analysis['distances_to_active_site'][active_res] = float('inf')
            
            # Pathways to active site
            pathways_to_active = self.find_allosteric_pathways([mut_site], self.active_site, interaction_type, threshold)
            site_analysis['pathways_to_active'] = pathways_to_active
            
            # Inter-chain communication (if applicable)
            if mut_site <= 99:  # Chain A
                partner_site = mut_site + 100  # Corresponding residue in chain B
                if partner_site in baseline_network:
                    try:
                        interchain_dist = nx.shortest_path_length(baseline_network, mut_site, partner_site, weight='weight')
                        site_analysis['interchain_distance'] = interchain_dist
                    except nx.NetworkXNoPath:
                        site_analysis['interchain_distance'] = float('inf')
            
            mutation_analysis[mut_site] = site_analysis
        
        return mutation_analysis

    def compare_wild_type_mutant_networks(self, wt_rin, mutant_rin, interaction_type='distance', threshold=0.2):
        """
        Compare network properties between wild-type and mutant systems
        
        Parameters:
        -----------
        wt_rin : HIVProteaseRIN
            Wild-type analysis object
        mutant_rin : HIVProteaseRIN  
            Mutant analysis object
        """
        print("Comparing WT and mutant networks...")
        
        wt_metrics = wt_rin.compute_hiv_specific_metrics(interaction_type, threshold)
        mut_metrics = mutant_rin.compute_hiv_specific_metrics(interaction_type, threshold)
        
        comparison = {}
        
        # Basic network property changes
        comparison['network_changes'] = {
            'density_change': mut_metrics['density'] - wt_metrics['density'],
            'clustering_change': mut_metrics['clustering_coefficient'] - wt_metrics['clustering_coefficient'],
            'path_length_change': mut_metrics['average_path_length'] - wt_metrics['average_path_length'],
            'modularity_change': mut_metrics['modularity'] - wt_metrics['modularity']
        }
        
        # Centrality changes for important residues
        important_residues = set(self.active_site + self.flap_regions['flap_tips'] + 
                               self.primary_resistance_sites + self.secondary_resistance_sites)
        
        comparison['centrality_changes'] = {}
        for res in important_residues:
            if res in wt_metrics['betweenness_centrality'] and res in mut_metrics['betweenness_centrality']:
                comparison['centrality_changes'][res] = {
                    'betweenness_change': mut_metrics['betweenness_centrality'][res] - wt_metrics['betweenness_centrality'][res],
                    'closeness_change': mut_metrics['closeness_centrality'][res] - wt_metrics['closeness_centrality'][res],
                    'eigenvector_change': mut_metrics['eigenvector_centrality'][res] - wt_metrics['eigenvector_centrality'][res]
                }
        
        # Pathway changes
        comparison['pathway_changes'] = {}
        wt_pathways = wt_metrics['flap_active_pathways']
        mut_pathways = mut_metrics['flap_active_pathways']
        
        for pathway_name in wt_pathways:
            if pathway_name in mut_pathways:
                wt_dist = wt_pathways[pathway_name]['shortest_distance']
                mut_dist = mut_pathways[pathway_name]['shortest_distance']
                comparison['pathway_changes'][pathway_name] = {
                    'distance_change': mut_dist - wt_dist,
                    'wt_paths': len(wt_pathways[pathway_name]['paths']),
                    'mut_paths': len(mut_pathways[pathway_name]['paths'])
                }
        
        # Community structure changes
        comparison['community_changes'] = {
            'wt_communities': len(wt_metrics['communities']),
            'mut_communities': len(mut_metrics['communities']),
            'modularity_change': mut_metrics['modularity'] - wt_metrics['modularity']
        }
        
        return comparison

    def create_hiv_specific_visualizations(self, prefix, metrics):
        """
        Create HIV protease-specific visualizations
        """
        # 1. Resistance site centrality plot
        self._plot_resistance_site_centrality(prefix, metrics)
        
        # 2. Flap-active site pathway visualization
        self._plot_flap_pathways(prefix, metrics)
        
        # 3. Inter-domain communication heatmap
        self._plot_interdomain_communication(prefix, metrics)
        
        # 4. Mutation effect radar plot
        self._plot_mutation_effects(prefix, metrics)

    def _plot_resistance_site_centrality(self, prefix, metrics):
        """Plot centrality of resistance sites"""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        centrality_types = ['betweenness_centrality', 'closeness_centrality', 
                           'eigenvector_centrality', 'degree_centrality']
        titles = ['Betweenness', 'Closeness', 'Eigenvector', 'Degree']
        
        for i, (cent_type, title) in enumerate(zip(centrality_types, titles)):
            ax = axes[i//2, i%2]
            
            # Primary resistance sites
            primary_vals = [metrics[cent_type].get(res, 0) for res in self.primary_resistance_sites]
            primary_labels = [str(res) for res in self.primary_resistance_sites]
            
            # Secondary resistance sites  
            secondary_vals = [metrics[cent_type].get(res, 0) for res in self.secondary_resistance_sites]
            secondary_labels = [str(res) for res in self.secondary_resistance_sites]
            
            x_prim = np.arange(len(primary_vals))
            x_sec = np.arange(len(secondary_vals)) + len(primary_vals) + 1
            
            ax.bar(x_prim, primary_vals, label='Primary', color='red', alpha=0.7)
            ax.bar(x_sec, secondary_vals, label='Secondary', color='blue', alpha=0.7)
            
            ax.set_title(f'{title} Centrality')
            ax.set_xlabel('Resistance Sites')
            ax.set_ylabel('Centrality Value')
            ax.set_xticks(np.concatenate([x_prim, x_sec]))
            ax.set_xticklabels(primary_labels + secondary_labels, rotation=45)
            ax.legend()
        
        plt.tight_layout()
        plt.savefig(f"{prefix}_resistance_centrality.png", dpi=300, bbox_inches='tight')
        plt.close()

    def _plot_flap_pathways(self, prefix, metrics):
        """Plot flap-to-active site pathway analysis"""
        pathways = metrics['flap_active_pathways']
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        pathway_names = []
        shortest_distances = []
        num_paths = []
        
        for pathway_name, pathway_data in pathways.items():
            pathway_names.append(pathway_name.replace('_to_', ' → '))
            shortest_distances.append(pathway_data['shortest_distance'])
            num_paths.append(len(pathway_data['paths']))
        
        x = np.arange(len(pathway_names))
        
        # Plot distances
        ax2 = ax.twinx()
        bars1 = ax.bar(x - 0.2, shortest_distances, 0.4, label='Shortest Distance', color='blue', alpha=0.7)
        bars2 = ax2.bar(x + 0.2, num_paths, 0.4, label='Number of Paths', color='red', alpha=0.7)
        
        ax.set_xlabel('Pathways')
        ax.set_ylabel('Shortest Distance (network units)', color='blue')
        ax2.set_ylabel('Number of Paths', color='red')
        
        ax.set_title('Flap-to-Active Site Communication Pathways')
        ax.set_xticks(x)
        ax.set_xticklabels(pathway_names, rotation=45, ha='right')
        
        # Combine legends
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left')
        
        plt.tight_layout()
        plt.savefig(f"{prefix}_flap_pathways.png", dpi=300, bbox_inches='tight')
        plt.close()

    def _plot_interdomain_communication(self, prefix, metrics):
        """Plot inter-domain communication heatmap"""
        interdomain_data = metrics['interdomain_communication']['interdomain_edges']
        
        # Create matrix
        domains = list(self.domains.keys())
        n_domains = len(domains)
        comm_matrix = np.zeros((n_domains, n_domains))
        
        for i, domain1 in enumerate(domains):
            for j, domain2 in enumerate(domains):
                if i < j:
                    key = f"{domain1}_{domain2}"
                    if key in interdomain_data:
                        weight = interdomain_data[key]['avg_weight']
                        comm_matrix[i, j] = weight
                        comm_matrix[j, i] = weight  # Symmetric
        
        plt.figure(figsize=(10, 8))
        sns.heatmap(comm_matrix, 
                   xticklabels=domains, 
                   yticklabels=domains,
                   annot=True, 
                   fmt='.3f',
                   cmap='viridis',
                   cbar_kws={'label': 'Average Edge Weight'})
        
        plt.title('Inter-Domain Communication Strength')
        plt.xlabel('Structural Domain')
        plt.ylabel('Structural Domain')
        plt.xticks(rotation=45)
        plt.yticks(rotation=0)
        plt.tight_layout()
        plt.savefig(f"{prefix}_interdomain_communication.png", dpi=300, bbox_inches='tight')
        plt.close()

    def _plot_mutation_effects(self, prefix, metrics):
        """Create a summary plot of overall network effects"""
        # This would be customized based on specific mutation analysis results
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Plot 1: Network density vs clustering
        ax = axes[0, 0]
        ax.scatter([metrics['density']], [metrics['clustering_coefficient']], 
                  s=100, c='red', alpha=0.7, label='Current System')
        ax.set_xlabel('Network Density')
        ax.set_ylabel('Clustering Coefficient')
        ax.set_title('Network Topology Properties')
        ax.legend()
        
        # Plot 2: Centrality distribution
        ax = axes[0, 1]
        centrality_vals = list(metrics['betweenness_centrality'].values())
        ax.hist(centrality_vals, bins=20, alpha=0.7, color='blue')
        ax.set_xlabel('Betweenness Centrality')
        ax.set_ylabel('Frequency')
        ax.set_title('Centrality Distribution')
        
        # Plot 3: Path lengths distribution
        ax = axes[1, 0]
        if 'shortest_path_lengths' in metrics:
            path_lengths = []
            for source_paths in metrics['shortest_path_lengths'].values():
                path_lengths.extend(source_paths.values())
            ax.hist(path_lengths, bins=20, alpha=0.7, color='green')
            ax.set_xlabel('Shortest Path Length')
            ax.set_ylabel('Frequency')
            ax.set_title('Path Length Distribution')
        
        # Plot 4: Community sizes
        ax = axes[1, 1]
        if metrics['communities']:
            community_sizes = [len(comm) for comm in metrics['communities']]
            ax.bar(range(len(community_sizes)), community_sizes, alpha=0.7, color='purple')
            ax.set_xlabel('Community ID')
            ax.set_ylabel('Community Size')
            ax.set_title('Community Structure')
        
        plt.tight_layout()
        plt.savefig(f"{prefix}_network_summary.png", dpi=300, bbox_inches='tight')
        plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="HIV Protease-specific residue interaction network analysis."
    )
    
    parser.add_argument("-t", "--topology", required=True,
                       help="Topology file")
    parser.add_argument("-x", "--trajectory", required=True,
                       help="Trajectory file")
    parser.add_argument("-o", "--output", default="hiv_protease_rin",
                       help="Output prefix")
    parser.add_argument("--threshold", type=float, default=0.2,
                       help="Contact frequency threshold")
    parser.add_argument("--mutation-sites", nargs='+', type=int,
                       help="Specific mutation sites to analyze")
    parser.add_argument("--compare-wt", 
                       help="Wild-type trajectory for comparison")
    parser.add_argument("--segments", type=int, default=5,
                       help="Number of trajectory segments")
    
    args = parser.parse_args()
    
    try:
        # Initialize HIV-specific RIN analysis
        hiv_rin = HIVProteaseRIN(args.topology, args.trajectory)
        
        # Compute contact maps
        print("Computing HIV protease contact maps...")
        cutoffs = {'all_atom': 4.5, 'ca_only': 8.0}
        contact_maps = hiv_rin.compute_multi_contact_maps(cutoffs, ['distance', 'hbond'])
        
        # HIV-specific network analysis
        print("Computing HIV-specific network metrics...")
        hiv_metrics = hiv_rin.compute_hiv_specific_metrics('distance', args.threshold)
        
        # Mutation analysis
        if args.mutation_sites:
            print("Analyzing mutation effects...")
            mutation_effects = hiv_rin.analyze_mutation_effects(args.mutation_sites)
            
            # Save mutation analysis
            import json
            with open(f"{args.output}_mutation_effects.json", 'w') as f:
                # Convert numpy types for JSON serialization
                def convert_numpy(obj):
                    if isinstance(obj, np.integer):
                        return int(obj)
                    elif isinstance(obj, np.floating):
                        return float(obj)
                    elif isinstance(obj, np.ndarray):
                        return obj.tolist()
                    return obj
                
                serializable_effects = {}
                for k, v in mutation_effects.items():
                    serializable_effects[k] = {
                        key: convert_numpy(val) for key, val in v.items() 
                        if not isinstance(val, dict) or key == 'distances_to_active_site'
                    }
                
                json.dump(serializable_effects, f, indent=2)
        
        # Wild-type comparison
        if args.compare_wt:
            print("Comparing with wild-type...")
            wt_rin = HIVProteaseRIN(args.topology, args.compare_wt)
            wt_contact_maps = wt_rin.compute_multi_contact_maps(cutoffs, ['distance'])
            
            comparison = hiv_rin.compare_wild_type_mutant_networks(wt_rin, hiv_rin)
            
            # Save comparison
            with open(f"{args.output}_wt_comparison.json", 'w') as f:
                json.dump(comparison, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
        
        # Segment analysis
        if args.segments > 1:
            print(f"Analyzing {args.segments} trajectory segments...")
            segment_results = hiv_rin.analyze_sub_trajectories(args.segments)
            
            import pickle
            with open(f"{args.output}_segment_analysis.pkl", 'wb') as f:
                pickle.dump(segment_results, f)
        
        # Save comprehensive outputs
        print("Saving outputs...")
        hiv_rin.save_comprehensive_outputs(args.output, hiv_metrics)
        hiv_rin.create_hiv_specific_visualizations(args.output, hiv_metrics)
        hiv_rin.create_pymol_visualization_script(args.output, hiv_metrics)
        
        # Print summary
        print(f"\nHIV Protease RIN Analysis Complete")
        print(f"Network: {hiv_metrics['n_nodes']} nodes, {hiv_metrics['n_edges']} edges")
        print(f"Density: {hiv_metrics['density']:.4f}")
        print(f"Modularity: {hiv_metrics['modularity']:.4f}")
        print(f"Flap-active pathways: {len(hiv_metrics['flap_active_pathways'])}")
        
        if args.mutation_sites:
            print(f"Analyzed {len(args.mutation_sites)} mutation sites")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
