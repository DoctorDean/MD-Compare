#!/usr/bin/env python3

"""
MD-Compare: A comprehensive toolkit for comparing molecular dynamics simulations

This module provides the core classes for loading, analyzing, and comparing
MD simulations using various network analysis methods.
"""

import os
import sys
import numpy as np
import pickle
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any
from dataclasses import dataclass
from abc import ABC, abstractmethod

try:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance, distance_array
    from MDAnalysis.analysis import contacts, align, rms
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
    import matplotlib.pyplot as plt
    import seaborn as sns
except ImportError:
    print("Warning: scipy, matplotlib, or seaborn not available. Some features will be limited.")


# =====================================================
# DATA STRUCTURES AND CONFIGURATION
# =====================================================

@dataclass
class SimulationConfig:
    """Configuration for a single MD simulation"""
    name: str
    topology: str
    trajectory: str
    selection: str = "protein and not name H*"
    description: str = ""
    metadata: Dict[str, Any] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}

@dataclass
class AnalysisConfig:
    """Configuration for network analysis parameters"""
    cutoffs: Dict[str, float] = None
    interaction_types: List[str] = None
    threshold: float = 0.2
    timeout_seconds: int = 300
    segments: int = 5
    preprocess: bool = True
    align_selection: str = "name CA"
    center_selection: str = "protein"
    
    def __post_init__(self):
        if self.cutoffs is None:
            self.cutoffs = {'all_atom': 4.5, 'ca_only': 8.0}
        if self.interaction_types is None:
            self.interaction_types = ['distance']

@dataclass
class NetworkMetrics:
    """Container for network analysis results"""
    n_nodes: int = 0
    n_edges: int = 0
    density: float = 0.0
    average_path_length: float = 0.0
    diameter: float = 0.0
    clustering_coefficient: float = 0.0
    modularity: float = 0.0
    communities: List[List[str]] = None
    betweenness_centrality: Dict[str, float] = None
    closeness_centrality: Dict[str, float] = None
    eigenvector_centrality: Dict[str, float] = None
    degree_centrality: Dict[str, float] = None
    network: nx.Graph = None
    
    def __post_init__(self):
        if self.communities is None:
            self.communities = []
        if self.betweenness_centrality is None:
            self.betweenness_centrality = {}
        if self.closeness_centrality is None:
            self.closeness_centrality = {}
        if self.eigenvector_centrality is None:
            self.eigenvector_centrality = {}
        if self.degree_centrality is None:
            self.degree_centrality = {}


# =====================================================
# CORE MD SIMULATION CLASS
# =====================================================

class MDSimulation:
    """
    Core class for handling individual MD simulation data
    """
    
    def __init__(self, config: SimulationConfig):
        """
        Initialize MD simulation from configuration
        
        Parameters:
        -----------
        config : SimulationConfig
            Configuration containing simulation details
        """
        self.config = config
        self.name = config.name
        self._universe = None
        self._atoms = None
        self._chain_info = {}
        self._residue_mapping = None
        
        # Analysis results storage
        self.contact_maps = {}
        self.network_metrics = None
        self.segment_results = {}
        
        print(f"Initializing MD simulation: {self.name}")
        
    def load(self) -> bool:
        """
        Load the MD simulation data
        
        Returns:
        --------
        bool : Success status
        """
        try:
            # Detect file format and load appropriately
            if self._is_desmond_system():
                print(f"Detected Desmond system for {self.name}")
                self._universe = mda.Universe(self.config.topology, self.config.trajectory)
            else:
                self._universe = mda.Universe(self.config.topology, self.config.trajectory)
                
            self._atoms = self._universe.select_atoms(self.config.selection)
            
            if len(self._atoms) == 0:
                raise RuntimeError(f"No atoms selected with '{self.config.selection}'")
            
            # Setup residue mapping
            self._setup_residue_mapping()
            
            print(f"Loaded {self.name}: {self.n_chains} chains, {self.n_residues} residues, {len(self._atoms)} atoms")
            print(f"Trajectory: {len(self._universe.trajectory)} frames")
            
            return True
            
        except Exception as e:
            print(f"Failed to load {self.name}: {e}")
            return False
    
    def _is_desmond_system(self) -> bool:
        """Check if this is a Desmond system"""
        desmond_exts = ['.cms', '.dms', '.dtr']
        return any(self.config.topology.endswith(ext) or 
                  self.config.trajectory.endswith(ext) for ext in desmond_exts)
    
    def _setup_residue_mapping(self):
        """Setup proper residue mapping for multi-chain systems"""
        # Get chain information
        chain_ids = self._atoms.chainIDs
        unique_chains = np.unique(chain_ids)
        
        self.unique_residue_keys = []
        self._chain_info = {}
        
        for chain in unique_chains:
            chain_atoms = self._atoms.select_atoms(f"chainid {chain}")
            chain_resids = np.unique(chain_atoms.resids)
            
            self._chain_info[chain] = {
                'residues': chain_resids,
                'n_residues': len(chain_resids),
                'atom_count': len(chain_atoms),
                'residue_range': (chain_resids[0], chain_resids[-1])
            }
            
            # Create unique keys for this chain
            for resid in chain_resids:
                unique_key = f"{chain}_{resid}"
                self.unique_residue_keys.append(unique_key)
        
        self.unique_residue_keys = np.array(self.unique_residue_keys)
        
        # Create mapping from atom index to unique residue index
        self._atom_to_residue_map = np.zeros(len(self._atoms), dtype=int)
        key_to_index = {key: i for i, key in enumerate(self.unique_residue_keys)}
        
        for i, atom in enumerate(self._atoms):
            unique_key = f"{atom.chainID}_{atom.resid}"
            if unique_key in key_to_index:
                self._atom_to_residue_map[i] = key_to_index[unique_key]
    
    @property
    def universe(self) -> mda.Universe:
        """Get the MDAnalysis Universe"""
        return self._universe
    
    @property
    def atoms(self) -> mda.AtomGroup:
        """Get the selected atoms"""
        return self._atoms
    
    @property
    def n_chains(self) -> int:
        """Number of chains in the system"""
        return len(self._chain_info)
    
    @property
    def n_residues(self) -> int:
        """Total number of residues"""
        return len(self.unique_residue_keys)
    
    @property
    def n_frames(self) -> int:
        """Number of trajectory frames"""
        return len(self._universe.trajectory) if self._universe else 0
    
    @property
    def chain_info(self) -> Dict:
        """Information about chains in the system"""
        return self._chain_info
    
    def get_system_summary(self) -> Dict[str, Any]:
        """Get a summary of the system properties"""
        return {
            'name': self.name,
            'n_chains': self.n_chains,
            'n_residues': self.n_residues,
            'n_atoms': len(self._atoms) if self._atoms else 0,
            'n_frames': self.n_frames,
            'chain_info': self.chain_info,
            'description': self.config.description,
            'metadata': self.config.metadata
        }


# =====================================================
# NETWORK ANALYSIS ENGINE
# =====================================================

class NetworkAnalyzer:
    """
    Engine for performing network analysis on MD simulations
    """
    
    def __init__(self, config: AnalysisConfig):
        """
        Initialize network analyzer
        
        Parameters:
        -----------
        config : AnalysisConfig
            Analysis configuration parameters
        """
        self.config = config
        self._preprocessor = None
    
    def compute_contact_maps(self, simulation: MDSimulation, 
                           start_frame: int = 0, end_frame: Optional[int] = None) -> Dict[str, np.ndarray]:
        """
        Compute contact maps for a simulation
        
        Parameters:
        -----------
        simulation : MDSimulation
            The simulation to analyze
        start_frame : int
            Starting frame for analysis
        end_frame : int, optional
            Ending frame (None = all frames)
            
        Returns:
        --------
        Dict[str, np.ndarray] : Contact maps by interaction type
        """
        print(f"Computing contact maps for {simulation.name}")
        
        contact_maps = {}
        for itype in self.config.interaction_types:
            contact_maps[itype] = np.zeros((simulation.n_residues, simulation.n_residues))
        
        if end_frame is None:
            end_frame = simulation.n_frames
        
        n_frames_analyzed = 0
        
        # Setup preprocessing if enabled
        if self.config.preprocess:
            self._setup_preprocessing(simulation)
        
        # Analyze frames
        for frame_idx in range(start_frame, end_frame):
            if frame_idx % 100 == 0:
                print(f"  Frame {frame_idx}/{end_frame} ({100*frame_idx/end_frame:.1f}%)")
            
            simulation.universe.trajectory[frame_idx]
            
            # Apply preprocessing
            if self.config.preprocess and self._preprocessor:
                self._preprocessor.center_system()
                self._preprocessor.align_to_reference()
            
            # Compute contacts
            if 'distance' in self.config.interaction_types:
                self._compute_distance_contacts(simulation, contact_maps['distance'])
            
            if 'hbond' in self.config.interaction_types:
                self._compute_hydrogen_bonds(simulation, contact_maps['hbond'])
                
            if 'salt_bridge' in self.config.interaction_types:
                self._compute_salt_bridges(simulation, contact_maps['salt_bridge'])
            
            n_frames_analyzed += 1
        
        # Normalize by number of frames
        for itype in contact_maps:
            contact_maps[itype] /= n_frames_analyzed
        
        # Store results
        simulation.contact_maps = contact_maps
        
        print(f"Contact maps computed over {n_frames_analyzed} frames")
        return contact_maps
    
    def _setup_preprocessing(self, simulation: MDSimulation):
        """Setup MD preprocessing for the simulation"""
        from .preprocessing import MDPreprocessor
        
        self._preprocessor = MDPreprocessor(
            simulation.universe,
            align_selection=self.config.align_selection,
            center_selection=self.config.center_selection
        )
        self._preprocessor.setup_reference(frame=0)
    
    def _compute_distance_contacts(self, simulation: MDSimulation, contact_matrix: np.ndarray):
        """Compute distance-based contacts"""
        coords = simulation.atoms.positions
        
        pairs = capped_distance(
            coords, coords,
            max_cutoff=self.config.cutoffs['all_atom'],
            return_distances=False
        )
        
        if len(pairs) > 0:
            atom_i, atom_j = pairs[:, 0], pairs[:, 1]
            res_i = simulation._atom_to_residue_map[atom_i]
            res_j = simulation._atom_to_residue_map[atom_j]
            
            # Only count inter-residue contacts
            mask = res_i != res_j
            res_i, res_j = res_i[mask], res_j[mask]
            
            contact_matrix[res_i, res_j] += 1
            contact_matrix[res_j, res_i] += 1
    
    def _compute_hydrogen_bonds(self, simulation: MDSimulation, contact_matrix: np.ndarray):
        """Compute hydrogen bond contacts"""
        # Implementation for H-bond detection
        pass
    
    def _compute_salt_bridges(self, simulation: MDSimulation, contact_matrix: np.ndarray):
        """Compute salt bridge contacts"""
        # Implementation for salt bridge detection
        pass
    
    def compute_network_metrics(self, simulation: MDSimulation, 
                               interaction_type: str = 'distance') -> NetworkMetrics:
        """
        Compute comprehensive network metrics
        
        Parameters:
        -----------
        simulation : MDSimulation
            The simulation with computed contact maps
        interaction_type : str
            Type of interaction network to analyze
            
        Returns:
        --------
        NetworkMetrics : Comprehensive network analysis results
        """
        print(f"Computing network metrics for {simulation.name}")
        
        if interaction_type not in simulation.contact_maps:
            raise ValueError(f"Contact map for '{interaction_type}' not found. Run compute_contact_maps first.")
        
        contact_matrix = simulation.contact_maps[interaction_type]
        
        # Build network
        G = self._build_network(simulation, contact_matrix)
        
        # Initialize metrics
        metrics = NetworkMetrics()
        metrics.network = G
        metrics.n_nodes = G.number_of_nodes()
        metrics.n_edges = G.number_of_edges()
        metrics.density = nx.density(G)
        
        if metrics.n_nodes == 0:
            print("Warning: Empty network")
            simulation.network_metrics = metrics
            return metrics
        
        # Compute centrality measures with timeout protection
        metrics = self._compute_centrality_measures(G, metrics)
        
        # Community detection
        metrics = self._detect_communities(G, metrics)
        
        # Path analysis
        metrics = self._analyze_paths(G, metrics)
        
        # Clustering
        try:
            metrics.clustering_coefficient = nx.average_clustering(G, weight='weight')
        except:
            metrics.clustering_coefficient = 0.0
        
        # Store results
        simulation.network_metrics = metrics
        
        print(f"Network analysis complete: {metrics.n_nodes} nodes, {metrics.n_edges} edges")
        return metrics
    
    def _build_network(self, simulation: MDSimulation, contact_matrix: np.ndarray) -> nx.Graph:
        """Build NetworkX graph from contact matrix"""
        G = nx.Graph()
        
        # Add nodes
        for i, res_key in enumerate(simulation.unique_residue_keys):
            chain, resid = res_key.split('_')
            G.add_node(res_key, chain=chain, resid=int(resid), index=i)
        
        # Add edges
        for i in range(simulation.n_residues):
            for j in range(i + 1, simulation.n_residues):
                weight = contact_matrix[i, j]
                if weight >= self.config.threshold:
                    G.add_edge(
                        simulation.unique_residue_keys[i],
                        simulation.unique_residue_keys[j],
                        weight=float(weight)
                    )
        
        return G
    
    def _compute_centrality_measures(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """Compute centrality measures with timeout protection"""
        from .utils import timeout_handler
        
        try:
            with timeout_handler(self.config.timeout_seconds):
                metrics.betweenness_centrality = nx.betweenness_centrality(
                    G, weight='weight', k=min(100, G.number_of_nodes())
                )
                metrics.closeness_centrality = nx.closeness_centrality(G, distance='weight')
                metrics.degree_centrality = nx.degree_centrality(G)
        except Exception as e:
            print(f"Centrality computation failed: {e}")
            metrics.betweenness_centrality = {node: 0 for node in G.nodes()}
            metrics.closeness_centrality = {node: 0 for node in G.nodes()}
            metrics.degree_centrality = {node: 0 for node in G.nodes()}
        
        # Eigenvector centrality (often problematic)
        try:
            with timeout_handler(60):
                metrics.eigenvector_centrality = nx.eigenvector_centrality(
                    G, weight='weight', max_iter=1000
                )
        except Exception:
            metrics.eigenvector_centrality = metrics.degree_centrality.copy()
        
        return metrics
    
    def _detect_communities(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """Detect network communities"""
        try:
            communities = nx.community.greedy_modularity_communities(G, weight='weight')
            metrics.communities = [list(c) for c in communities]
            metrics.modularity = nx.community.modularity(G, communities, weight='weight')
        except Exception:
            metrics.communities = []
            metrics.modularity = 0.0
        
        return metrics
    
    def _analyze_paths(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """Analyze network paths"""
        try:
            if nx.is_connected(G):
                if G.number_of_nodes() > 500:
                    # Sample-based approximation for large networks
                    sample_nodes = np.random.choice(list(G.nodes()), min(50, len(G.nodes())), replace=False)
                    sample_lengths = []
                    for node in sample_nodes:
                        lengths = nx.single_source_shortest_path_length(G, node, weight='weight')
                        sample_lengths.extend([l for l in lengths.values() if l > 0])
                    
                    metrics.average_path_length = np.mean(sample_lengths) if sample_lengths else 0
                    metrics.diameter = max(sample_lengths) if sample_lengths else 0
                else:
                    metrics.average_path_length = nx.average_shortest_path_length(G, weight='weight')
                    metrics.diameter = nx.diameter(G)
            else:
                # Disconnected network - analyze largest component
                components = list(nx.connected_components(G))
                if components:
                    largest_cc = max(components, key=len)
                    G_cc = G.subgraph(largest_cc)
                    if len(G_cc) > 1:
                        metrics.average_path_length = nx.average_shortest_path_length(G_cc, weight='weight')
                        metrics.diameter = nx.diameter(G_cc)
        except Exception:
            metrics.average_path_length = float('inf')
            metrics.diameter = float('inf')
        
        return metrics
    
    def analyze_trajectory_segments(self, simulation: MDSimulation) -> Dict[str, Any]:
        """
        Analyze trajectory in segments for statistical robustness
        
        Parameters:
        -----------
        simulation : MDSimulation
            The simulation to analyze
            
        Returns:
        --------
        Dict : Results from segmented analysis
        """
        print(f"Analyzing {self.config.segments} trajectory segments for {simulation.name}")
        
        segment_size = simulation.n_frames // self.config.segments
        segment_results = {}
        
        for seg in range(self.config.segments):
            print(f"Segment {seg + 1}/{self.config.segments}")
            
            start_frame = seg * segment_size
            end_frame = start_frame + segment_size
            
            # Compute contact maps for this segment
            segment_contacts = self.compute_contact_maps(simulation, start_frame, end_frame)
            
            # Compute network metrics for this segment
            segment_metrics = self.compute_network_metrics(simulation, 'distance')
            
            segment_results[f'segment_{seg}'] = {
                'frame_range': (start_frame, end_frame),
                'contact_maps': segment_contacts,
                'network_metrics': segment_metrics
            }
        
        simulation.segment_results = segment_results
        return segment_results


# =====================================================
# COMPARISON ENGINE
# =====================================================

class MDComparator:
    """
    Engine for comparing multiple MD simulations
    """
    
    def __init__(self, simulations: List[MDSimulation], analyzer: NetworkAnalyzer):
        """
        Initialize MD comparator
        
        Parameters:
        -----------
        simulations : List[MDSimulation]
            List of simulations to compare
        analyzer : NetworkAnalyzer
            Network analyzer instance
        """
        self.simulations = {sim.name: sim for sim in simulations}
        self.analyzer = analyzer
        self.comparison_results = {}
    
    def add_simulation(self, simulation: MDSimulation):
        """Add a simulation to the comparison set"""
        self.simulations[simulation.name] = simulation
    
    def remove_simulation(self, name: str):
        """Remove a simulation from the comparison set"""
        if name in self.simulations:
            del self.simulations[name]
    
    def compare_network_properties(self) -> Dict[str, Any]:
        """
        Compare basic network properties across simulations
        
        Returns:
        --------
        Dict : Comparison of network properties
        """
        comparison = {}
        
        for name, sim in self.simulations.items():
            if sim.network_metrics is None:
                print(f"Warning: No network metrics for {name}. Run analysis first.")
                continue
            
            comparison[name] = {
                'n_nodes': sim.network_metrics.n_nodes,
                'n_edges': sim.network_metrics.n_edges,
                'density': sim.network_metrics.density,
                'average_path_length': sim.network_metrics.average_path_length,
                'diameter': sim.network_metrics.diameter,
                'clustering_coefficient': sim.network_metrics.clustering_coefficient,
                'modularity': sim.network_metrics.modularity,
                'n_communities': len(sim.network_metrics.communities)
            }
        
        self.comparison_results['network_properties'] = comparison
        return comparison
    
    def compare_centrality_measures(self, residue_keys: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Compare centrality measures for specific residues
        
        Parameters:
        -----------
        residue_keys : List[str], optional
            Specific residues to compare (None = all common residues)
            
        Returns:
        --------
        Dict : Centrality comparison data
        """
        if residue_keys is None:
            # Find common residues across all simulations
            all_residues = [set(sim.unique_residue_keys) for sim in self.simulations.values()]
            common_residues = set.intersection(*all_residues) if all_residues else set()
            residue_keys = list(common_residues)
        
        centrality_comparison = {}
        centrality_types = ['betweenness_centrality', 'closeness_centrality', 
                           'eigenvector_centrality', 'degree_centrality']
        
        for cent_type in centrality_types:
            centrality_comparison[cent_type] = {}
            
            for res_key in residue_keys:
                centrality_comparison[cent_type][res_key] = {}
                
                for name, sim in self.simulations.items():
                    if sim.network_metrics and hasattr(sim.network_metrics, cent_type):
                        centralities = getattr(sim.network_metrics, cent_type)
                        centrality_comparison[cent_type][res_key][name] = centralities.get(res_key, 0)
        
        self.comparison_results['centrality_measures'] = centrality_comparison
        return centrality_comparison
    
    def compare_contact_patterns(self, interaction_type: str = 'distance') -> Dict[str, Any]:
        """
        Compare contact patterns between simulations
        
        Parameters:
        -----------
        interaction_type : str
            Type of contacts to compare
            
        Returns:
        --------
        Dict : Contact pattern comparison
        """
        contact_comparison = {}
        
        # Find common residue pairs
        common_residues = None
        for sim in self.simulations.values():
            if common_residues is None:
                common_residues = set(sim.unique_residue_keys)
            else:
                common_residues = common_residues.intersection(set(sim.unique_residue_keys))
        
        if not common_residues:
            print("Warning: No common residues found across simulations")
            return {}
        
        common_residues = list(common_residues)
        
        # Compare contact frequencies
        for name, sim in self.simulations.items():
            if interaction_type not in sim.contact_maps:
                continue
                
            contact_matrix = sim.contact_maps[interaction_type]
            
            # Extract relevant subset
            residue_indices = [list(sim.unique_residue_keys).index(res) for res in common_residues]
            subset_matrix = contact_matrix[np.ix_(residue_indices, residue_indices)]
            
            contact_comparison[name] = {
                'contact_matrix': subset_matrix,
                'total_contacts': np.sum(subset_matrix > 0) // 2,  # Symmetric matrix
                'average_contact_frequency': np.mean(subset_matrix[subset_matrix > 0]),
                'max_contact_frequency': np.max(subset_matrix)
            }
        
        self.comparison_results['contact_patterns'] = contact_comparison
        return contact_comparison
    
    def find_differential_contacts(self, sim1_name: str, sim2_name: str, 
                                  threshold_diff: float = 0.1) -> Dict[str, Any]:
        """
        Find contacts that differ significantly between two simulations
        
        Parameters:
        -----------
        sim1_name : str
            Name of first simulation
        sim2_name : str
            Name of second simulation
        threshold_diff : float
            Minimum difference in contact frequency to consider significant
            
        Returns:
        --------
        Dict : Differential contacts analysis
        """
        if sim1_name not in self.simulations or sim2_name not in self.simulations:
            raise ValueError(f"Simulation(s) not found")
        
        sim1 = self.simulations[sim1_name]
        sim2 = self.simulations[sim2_name]
        
        # Find common residues
        common_residues = set(sim1.unique_residue_keys).intersection(set(sim2.unique_residue_keys))
        common_residues = list(common_residues)
        
        # Get contact matrices
        contact1 = sim1.contact_maps.get('distance')
        contact2 = sim2.contact_maps.get('distance')
        
        if contact1 is None or contact2 is None:
            raise ValueError("Contact maps not available for both simulations")
        
        # Extract common residue subsets
        indices1 = [list(sim1.unique_residue_keys).index(res) for res in common_residues]
        indices2 = [list(sim2.unique_residue_keys).index(res) for res in common_residues]
        
        matrix1 = contact1[np.ix_(indices1, indices1)]
        matrix2 = contact2[np.ix_(indices2, indices2)]
        
        # Calculate differences
        diff_matrix = matrix2 - matrix1
        
        # Find significant differences
        significant_increases = []
        significant_decreases = []
        
        for i in range(len(common_residues)):
            for j in range(i + 1, len(common_residues)):
                diff = diff_matrix[i, j]
                if abs(diff) >= threshold_diff:
                    residue_pair = (common_residues[i], common_residues[j])
                    contact_data = {
                        'residue_pair': residue_pair,
                        'sim1_frequency': matrix1[i, j],
                        'sim2_frequency': matrix2[i, j],
                        'difference': diff
                    }
                    
                    if diff > 0:
                        significant_increases.append(contact_data)
                    else:
                        significant_decreases.append(contact_data)
        
        differential_results = {
            'comparison': f"{sim1_name} vs {sim2_name}",
            'threshold_difference': threshold_diff,
            'increased_contacts': sorted(significant_increases, key=lambda x: x['difference'], reverse=True),
            'decreased_contacts': sorted(significant_decreases, key=lambda x: x['difference']),
            'n_increased': len(significant_increases),
            'n_decreased': len(significant_decreases),
            'difference_matrix': diff_matrix,
            'common_residues': common_residues
        }
        
        return differential_results
    
    def generate_comparison_report(self) -> Dict[str, Any]:
        """
        Generate a comprehensive comparison report
        
        Returns:
        --------
        Dict : Complete comparison analysis
        """
        report = {
            'simulations': {name: sim.get_system_summary() for name, sim in self.simulations.items()},
            'analysis_config': {
                'cutoffs': self.analyzer.config.cutoffs,
                'threshold': self.analyzer.config.threshold,
                'interaction_types': self.analyzer.config.interaction_types,
                'segments': self.analyzer.config.segments
            }
        }
        
        # Add comparison results
        report.update(self.comparison_results)
        
        return report


# =====================================================
# OUTPUT MANAGER
# =====================================================

class OutputManager:
    """
    Manages output generation and file saving for MD comparisons
    """
    
    def __init__(self, output_dir: str = "md_compare_results"):
        """
        Initialize output manager
        
        Parameters:
        -----------
        output_dir : str
            Directory for saving outputs
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
    def save_simulation_results(self, simulation: MDSimulation, prefix: str = None):
        """Save results for a single simulation"""
        if prefix is None:
            prefix = simulation.name
        
        # Save contact maps
        for itype, contact_matrix in simulation.contact_maps.items():
            np.save(self.output_dir / f"{prefix}_{itype}_contacts.npy", contact_matrix)
            
            # Save as CSV if not too large
            if contact_matrix.size < 1000000:
                np.savetxt(
                    self.output_dir / f"{prefix}_{itype}_contacts.csv",
                    contact_matrix,
                    delimiter=",",
                    header=",".join(simulation.unique_residue_keys),
                    comments=""
                )
        
        # Save network
        if simulation.network_metrics and simulation.network_metrics.network:
            nx.write_graphml(
                simulation.network_metrics.network,
                self.output_dir / f"{prefix}_network.graphml"
            )
        
        # Save centrality data
        if simulation.network_metrics:
            self._save_centrality_data(simulation, prefix)
        
        # Save system info
        self._save_system_info(simulation, prefix)
    
    def save_comparison_results(self, comparator: MDComparator, prefix: str = "comparison"):
        """Save comparison results"""
        report = comparator.generate_comparison_report()
        
        # Save as JSON
        with open(self.output_dir / f"{prefix}_report.json", 'w') as f:
            json.dump(report, f, indent=2, default=str)
        
        # Save as pickle for full data
        with open(self.output_dir / f"{prefix}_full_results.pkl", 'wb') as f:
            pickle.dump(report, f)
    
    def _save_centrality_data(self, simulation: MDSimulation, prefix: str):
        """Save centrality measures as CSV"""
        if not simulation.network_metrics:
            return
            
        centrality_data = []
        for res_key in simulation.unique_residue_keys:
            chain, resid = res_key.split('_')
            centrality_data.append([
                res_key,
                chain,
                resid,
                simulation.network_metrics.betweenness_centrality.get(res_key, 0),
                simulation.network_metrics.closeness_centrality.get(res_key, 0),
                simulation.network_metrics.eigenvector_centrality.get(res_key, 0),
                simulation.network_metrics.degree_centrality.get(res_key, 0)
            ])
        
        np.savetxt(
            self.output_dir / f"{prefix}_centrality.csv",
            centrality_data,
            delimiter=",",
            header="residue_key,chain,resid,betweenness,closeness,eigenvector,degree",
            fmt='%s',
            comments=""
        )
    
    def _save_system_info(self, simulation: MDSimulation, prefix: str):
        """Save system information"""
        with open(self.output_dir / f"{prefix}_system_info.txt", 'w') as f:
            summary = simulation.get_system_summary()
            f.write(f"MD Simulation Analysis: {simulation.name}\n")
            f.write("=" * 50 + "\n\n")
            
            for key, value in summary.items():
                if key != 'chain_info':
                    f.write(f"{key}: {value}\n")
            
            f.write("\nChain Information:\n")
            for chain, info in summary['chain_info'].items():
                f.write(f"  Chain {chain}:\n")
                f.write(f"    Residues: {info['n_residues']}\n")
                f.write(f"    Range: {info['residue_range'][0]}-{info['residue_range'][1]}\n")
                f.write(f"    Atoms: {info['atom_count']}\n\n")


# =====================================================
# MAIN WORKFLOW CLASS
# =====================================================

class MDCompare:
    """
    Main workflow class for MD comparison analysis
    """
    
    def __init__(self, analysis_config: AnalysisConfig = None, output_dir: str = "md_compare_results"):
        """
        Initialize MD-Compare workflow
        
        Parameters:
        -----------
        analysis_config : AnalysisConfig, optional
            Configuration for analysis parameters
        output_dir : str
            Directory for saving results
        """
        self.analysis_config = analysis_config or AnalysisConfig()
        self.analyzer = NetworkAnalyzer(self.analysis_config)
        self.output_manager = OutputManager(output_dir)
        self.simulations = {}
        self.comparator = None
        
        print(f"MD-Compare initialized with output directory: {output_dir}")
    
    def add_simulation(self, config: SimulationConfig) -> bool:
        """
        Add a simulation to the comparison set
        
        Parameters:
        -----------
        config : SimulationConfig
            Simulation configuration
            
        Returns:
        --------
        bool : Success status
        """
        simulation = MDSimulation(config)
        success = simulation.load()
        
        if success:
            self.simulations[config.name] = simulation
            print(f"Added simulation: {config.name}")
        else:
            print(f"Failed to add simulation: {config.name}")
        
        return success
    
    def run_analysis(self, simulation_names: List[str] = None) -> Dict[str, Any]:
        """
        Run complete analysis on specified simulations
        
        Parameters:
        -----------
        simulation_names : List[str], optional
            Names of simulations to analyze (None = all)
            
        Returns:
        --------
        Dict : Analysis results summary
        """
        if simulation_names is None:
            simulation_names = list(self.simulations.keys())
        
        results_summary = {}
        
        for name in simulation_names:
            if name not in self.simulations:
                print(f"Warning: Simulation '{name}' not found")
                continue
            
            print(f"\nAnalyzing {name}...")
            simulation = self.simulations[name]
            
            # Compute contact maps
            contact_maps = self.analyzer.compute_contact_maps(simulation)
            
            # Compute network metrics
            network_metrics = self.analyzer.compute_network_metrics(simulation)
            
            # Trajectory segmentation analysis
            if self.analysis_config.segments > 1:
                segment_results = self.analyzer.analyze_trajectory_segments(simulation)
            
            # Save individual results
            self.output_manager.save_simulation_results(simulation)
            
            results_summary[name] = {
                'contact_maps_computed': list(contact_maps.keys()),
                'network_nodes': network_metrics.n_nodes,
                'network_edges': network_metrics.n_edges,
                'network_density': network_metrics.density
            }
        
        return results_summary
    
    def run_comparison(self, simulation_names: List[str] = None) -> Dict[str, Any]:
        """
        Run comparative analysis between simulations
        
        Parameters:
        -----------
        simulation_names : List[str], optional
            Names of simulations to compare (None = all)
            
        Returns:
        --------
        Dict : Comparison results
        """
        if simulation_names is None:
            simulation_names = list(self.simulations.keys())
        
        # Filter available simulations
        available_sims = [self.simulations[name] for name in simulation_names 
                         if name in self.simulations]
        
        if len(available_sims) < 2:
            raise ValueError("Need at least 2 simulations for comparison")
        
        # Create comparator
        self.comparator = MDComparator(available_sims, self.analyzer)
        
        # Run comparisons
        print("Running comparative analysis...")
        
        network_comparison = self.comparator.compare_network_properties()
        centrality_comparison = self.comparator.compare_centrality_measures()
        contact_comparison = self.comparator.compare_contact_patterns()
        
        # Save comparison results
        self.output_manager.save_comparison_results(self.comparator)
        
        return {
            'network_properties': network_comparison,
            'centrality_measures': centrality_comparison,
            'contact_patterns': contact_comparison
        }
    
    def run_full_workflow(self, configs: List[SimulationConfig]) -> Dict[str, Any]:
        """
        Run complete MD-Compare workflow
        
        Parameters:
        -----------
        configs : List[SimulationConfig]
            List of simulation configurations
            
        Returns:
        --------
        Dict : Complete analysis results
        """
        print("Starting MD-Compare full workflow...")
        
        # Add all simulations
        successful_simulations = []
        for config in configs:
            if self.add_simulation(config):
                successful_simulations.append(config.name)
        
        if len(successful_simulations) == 0:
            raise RuntimeError("No simulations loaded successfully")
        
        # Run individual analyses
        analysis_results = self.run_analysis(successful_simulations)
        
        # Run comparative analysis if multiple simulations
        comparison_results = {}
        if len(successful_simulations) > 1:
            comparison_results = self.run_comparison(successful_simulations)
        
        # Generate final report
        final_report = {
            'workflow_config': {
                'analysis_config': self.analysis_config.__dict__,
                'output_directory': str(self.output_manager.output_dir)
            },
            'simulations_analyzed': successful_simulations,
            'individual_analysis': analysis_results,
            'comparative_analysis': comparison_results
        }
        
        # Save final report
        with open(self.output_manager.output_dir / "workflow_summary.json", 'w') as f:
            json.dump(final_report, f, indent=2, default=str)
        
        print(f"\nMD-Compare workflow complete!")
        print(f"Results saved to: {self.output_manager.output_dir}")
        
        return final_report
