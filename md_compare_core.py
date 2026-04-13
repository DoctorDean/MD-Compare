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
    from scipy.linalg import svd
    from scipy.ndimage import gaussian_filter
    import matplotlib.pyplot as plt
    import seaborn as sns
    # Advanced network analysis
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path, connected_components
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
    from sklearn.cluster import SpectralClustering
    from collections import deque, defaultdict
except ImportError:
    print("Warning: scipy, matplotlib, seaborn, or sklearn not available. Some features will be limited.")


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
    # Dynamic analysis options
    compute_dccm: bool = True
    compute_pca: bool = True
    pca_components: int = 10
    dccm_selection: str = "name CA"
    # Energy landscape analysis options
    compute_energy_landscape: bool = True
    landscape_temperature: float = 310.0  # Kelvin
    landscape_bins: int = 50
    landscape_sigma: float = 1.0  # Gaussian smoothing
    landscape_epsilon: float = 1e-10  # Zero bin handling
    # Advanced network analysis options
    compute_communities: bool = True
    community_method: str = "leiden"  # leiden, louvain, spectral, hierarchical
    compute_paths: bool = True
    path_analysis_nodes: List[str] = None  # Specific nodes for detailed path analysis
    allosteric_analysis: bool = True
    allosteric_source_nodes: List[str] = None  # Source residues for allosteric analysis
    allosteric_target_nodes: List[str] = None  # Target residues for allosteric analysis
    
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
    # Dynamic analysis fields
    dccm_matrix: Optional[np.ndarray] = None
    pca_eigenvalues: Optional[np.ndarray] = None
    pca_eigenvectors: Optional[np.ndarray] = None
    pca_variance_explained: Optional[np.ndarray] = None
    principal_components: Optional[np.ndarray] = None
    # Energy landscape fields
    energy_landscape: Optional[np.ndarray] = None
    landscape_pc1_bins: Optional[np.ndarray] = None
    landscape_pc2_bins: Optional[np.ndarray] = None
    landscape_gradient_x: Optional[np.ndarray] = None
    landscape_gradient_y: Optional[np.ndarray] = None
    landscape_laplacian: Optional[np.ndarray] = None
    landscape_minima: Optional[List[Tuple[float, float, float]]] = None
    landscape_barriers: Optional[List[Dict[str, Any]]] = None
    # Advanced network analysis fields
    communities_detailed: Optional[Dict[str, Any]] = None
    community_modularity: Optional[float] = None
    community_node_assignments: Optional[Dict[str, int]] = None
    path_metrics: Optional[Dict[str, Any]] = None
    allosteric_pathways: Optional[Dict[str, Any]] = None
    network_robustness: Optional[Dict[str, float]] = None
    centrality_z_scores: Optional[Dict[str, Dict[str, float]]] = None
    
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
        from utils import MDPreprocessor
        
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
        
        # Advanced network analysis
        if self.config.compute_communities:
            metrics = self._advanced_community_detection(G, metrics)
        
        # Path analysis
        if self.config.compute_paths:
            metrics = self._compute_path_metrics(G, metrics)
        
        # Allosteric pathway analysis
        if self.config.allosteric_analysis:
            metrics = self._analyze_allosteric_pathways(G, metrics, simulation)
        
        # Network robustness analysis
        metrics = self._analyze_network_robustness(G, metrics)
        
        # Compute centrality z-scores for significance testing
        metrics = self._compute_centrality_significance(G, metrics)
        
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
    
    def compute_dynamic_analysis(self, simulation: MDSimulation) -> Dict[str, Any]:
        """
        Compute dynamic cross-correlation matrix (DCCM) and principal component analysis (PCA)
        
        Parameters:
        -----------
        simulation : MDSimulation
            The simulation to analyze
            
        Returns:
        --------
        Dict : Dynamic analysis results including DCCM and PCA
        """
        dynamic_results = {}
        
        if not self.config.compute_dccm and not self.config.compute_pca:
            return dynamic_results
            
        print("Computing dynamic analysis (DCCM and PCA)...")
        
        # Get atoms for analysis (typically CA atoms)
        atoms = simulation.universe.select_atoms(self.config.dccm_selection)
        
        if len(atoms) == 0:
            print(f"Warning: No atoms selected with '{self.config.dccm_selection}'")
            return dynamic_results
            
        print(f"Dynamic analysis on {len(atoms)} atoms")
        
        # Collect coordinates across trajectory
        coordinates = []
        n_frames = len(simulation.universe.trajectory)
        
        print("Collecting coordinates...")
        for i, ts in enumerate(simulation.universe.trajectory):
            if i % 100 == 0:
                print(f"  Frame {i}/{n_frames} ({100*i/n_frames:.1f}%)")
                
            # Apply preprocessing if enabled
            if self.config.preprocess and self._preprocessor:
                self._preprocessor.center_system()
                self._preprocessor.align_to_reference()
                
            coordinates.append(atoms.positions.copy())
        
        coordinates = np.array(coordinates)  # Shape: (n_frames, n_atoms, 3)
        
        # Compute DCCM
        if self.config.compute_dccm:
            print("Computing DCCM...")
            dccm_matrix = self._compute_dccm(coordinates)
            dynamic_results['dccm_matrix'] = dccm_matrix
            dynamic_results['dccm_atoms'] = len(atoms)
            
        # Compute PCA
        if self.config.compute_pca:
            print("Computing PCA...")
            pca_results = self._compute_pca(coordinates, self.config.pca_components)
            dynamic_results.update(pca_results)
            
            # Compute energy landscape from PC projections
            if self.config.compute_energy_landscape and 'principal_components' in pca_results:
                print("Computing energy landscape...")
                landscape_results = self._compute_energy_landscape(pca_results['principal_components'])
                dynamic_results.update(landscape_results)
            
        print("Dynamic analysis completed!")
        return dynamic_results
    
    def _compute_dccm(self, coordinates: np.ndarray) -> np.ndarray:
        """
        Compute Dynamic Cross-Correlation Matrix
        
        Parameters:
        -----------
        coordinates : np.ndarray
            Coordinates array (n_frames, n_atoms, 3)
            
        Returns:
        --------
        np.ndarray : DCCM matrix (n_atoms, n_atoms)
        """
        n_frames, n_atoms, _ = coordinates.shape
        
        # Calculate displacement vectors (coordinates - average)
        avg_coords = np.mean(coordinates, axis=0)
        displacements = coordinates - avg_coords[np.newaxis, :, :]
        
        # Initialize DCCM matrix
        dccm = np.zeros((n_atoms, n_atoms))
        
        # Compute cross-correlations
        for i in range(n_atoms):
            for j in range(i, n_atoms):
                # Displacement vectors for atoms i and j
                disp_i = displacements[:, i, :]  # Shape: (n_frames, 3)
                disp_j = displacements[:, j, :]  # Shape: (n_frames, 3)
                
                # Compute cross-correlation
                numerator = np.mean(np.sum(disp_i * disp_j, axis=1))
                
                # Compute normalization factors
                variance_i = np.mean(np.sum(disp_i**2, axis=1))
                variance_j = np.mean(np.sum(disp_j**2, axis=1))
                
                if variance_i > 0 and variance_j > 0:
                    correlation = numerator / np.sqrt(variance_i * variance_j)
                else:
                    correlation = 0.0
                    
                dccm[i, j] = correlation
                dccm[j, i] = correlation  # Symmetric matrix
                
        return dccm
    
    def _compute_pca(self, coordinates: np.ndarray, n_components: int) -> Dict[str, Any]:
        """
        Compute Principal Component Analysis
        
        Parameters:
        -----------
        coordinates : np.ndarray
            Coordinates array (n_frames, n_atoms, 3)
        n_components : int
            Number of principal components to compute
            
        Returns:
        --------
        Dict : PCA results including eigenvalues, eigenvectors, and variance explained
        """
        n_frames, n_atoms, _ = coordinates.shape
        
        # Reshape coordinates to (n_frames, n_atoms*3)
        coords_flat = coordinates.reshape(n_frames, -1)
        
        # Center the data
        mean_coords = np.mean(coords_flat, axis=0)
        centered_coords = coords_flat - mean_coords
        
        # Compute covariance matrix
        cov_matrix = np.cov(centered_coords.T)
        
        # Perform SVD
        try:
            eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
            
            # Sort by eigenvalue magnitude (descending)
            idx = np.argsort(eigenvalues)[::-1]
            eigenvalues = eigenvalues[idx]
            eigenvectors = eigenvectors[:, idx]
            
            # Take only requested number of components
            n_components = min(n_components, len(eigenvalues))
            eigenvalues = eigenvalues[:n_components]
            eigenvectors = eigenvectors[:, :n_components]
            
            # Compute variance explained
            total_variance = np.sum(eigenvalues)
            variance_explained = eigenvalues / total_variance if total_variance > 0 else np.zeros_like(eigenvalues)
            cumulative_variance = np.cumsum(variance_explained)
            
            # Project data onto principal components
            principal_components = np.dot(centered_coords, eigenvectors)
            
            pca_results = {
                'pca_eigenvalues': eigenvalues,
                'pca_eigenvectors': eigenvectors,
                'pca_variance_explained': variance_explained,
                'pca_cumulative_variance': cumulative_variance,
                'principal_components': principal_components,
                'pca_n_components': n_components,
                'pca_total_variance': total_variance
            }
            
            print(f"PCA: {n_components} components explain {cumulative_variance[-1]*100:.1f}% of variance")
            
        except Exception as e:
            print(f"PCA computation failed: {e}")
            pca_results = {
                'pca_eigenvalues': None,
                'pca_eigenvectors': None,
                'pca_variance_explained': None,
                'principal_components': None
            }
            
        return pca_results
    
    def _compute_energy_landscape(self, principal_components: np.ndarray) -> Dict[str, Any]:
        """
        Compute free energy landscape from PC1 and PC2 projections
        
        Uses the Boltzmann relation: G = -kB*T*ln(P) where P is probability density
        
        Parameters:
        -----------
        principal_components : np.ndarray
            PC projections (n_frames, n_components)
            
        Returns:
        --------
        Dict : Energy landscape analysis results
        """
        if principal_components.shape[1] < 2:
            print("Warning: Need at least 2 PCs for energy landscape analysis")
            return {}
            
        try:
            # Extract PC1 and PC2
            pc1 = principal_components[:, 0]
            pc2 = principal_components[:, 1]
            
            # Physical constants
            kB = 8.314462618e-3  # kJ/(mol·K) - Boltzmann constant
            T = self.config.landscape_temperature  # Temperature in Kelvin
            epsilon = self.config.landscape_epsilon
            
            # Create bins for histogram
            n_bins = self.config.landscape_bins
            pc1_bins = np.linspace(np.min(pc1), np.max(pc1), n_bins)
            pc2_bins = np.linspace(np.min(pc2), np.max(pc2), n_bins)
            
            print(f"Creating {n_bins}x{n_bins} energy landscape")
            print(f"PC1 range: {np.min(pc1):.2f} to {np.max(pc1):.2f}")
            print(f"PC2 range: {np.min(pc2):.2f} to {np.max(pc2):.2f}")
            
            # Compute 2D histogram
            H, pc1_edges, pc2_edges = np.histogram2d(pc1, pc2, bins=[pc1_bins, pc2_bins])
            H = H.T  # Transpose for correct orientation
            
            # Convert to probability density
            P = H / np.sum(H)
            
            # Handle zero bins (add small epsilon)
            P[P == 0] = epsilon
            
            # Compute free energy using Boltzmann relation
            G = -kB * T * np.log(P)
            
            # Set minimum energy to zero (relative energies)
            G -= np.min(G)
            
            # Apply Gaussian smoothing
            if self.config.landscape_sigma > 0:
                G = gaussian_filter(G, sigma=self.config.landscape_sigma)
                # Re-normalize after smoothing
                G -= np.min(G)
            
            print(f"Energy landscape: {np.min(G):.1f} to {np.max(G):.1f} kJ/mol")
            
            # Compute gradients (forces)
            gradient_y, gradient_x = np.gradient(G)
            
            # Compute Laplacian (curvature)
            laplacian = self._compute_laplacian(G)
            
            # Find local minima and barriers
            minima = self._find_energy_minima(G, pc1_edges, pc2_edges)
            barriers = self._analyze_energy_barriers(G, minima, pc1_edges, pc2_edges)
            
            landscape_results = {
                'energy_landscape': G,
                'landscape_pc1_bins': pc1_edges,
                'landscape_pc2_bins': pc2_edges,
                'landscape_gradient_x': gradient_x,
                'landscape_gradient_y': gradient_y,
                'landscape_laplacian': laplacian,
                'landscape_minima': minima,
                'landscape_barriers': barriers,
                'landscape_temperature': T,
                'landscape_energy_range': [np.min(G), np.max(G)],
                'landscape_total_frames': len(pc1)
            }
            
            print(f"Found {len(minima)} energy minima")
            if barriers:
                print(f"Identified {len(barriers)} energy barriers")
            
            return landscape_results
            
        except Exception as e:
            print(f"Energy landscape computation failed: {e}")
            return {}
    
    def _compute_laplacian(self, energy_surface: np.ndarray) -> np.ndarray:
        """
        Compute Laplacian (second derivative) of energy surface
        
        Parameters:
        -----------
        energy_surface : np.ndarray
            2D energy landscape
            
        Returns:
        --------
        np.ndarray : Laplacian matrix
        """
        # Compute second derivatives using finite differences
        # Laplacian = d²G/dx² + d²G/dy²
        
        # Second derivative in x direction
        d2dx2 = np.zeros_like(energy_surface)
        d2dx2[:, 1:-1] = energy_surface[:, :-2] - 2*energy_surface[:, 1:-1] + energy_surface[:, 2:]
        
        # Second derivative in y direction  
        d2dy2 = np.zeros_like(energy_surface)
        d2dy2[1:-1, :] = energy_surface[:-2, :] - 2*energy_surface[1:-1, :] + energy_surface[2:, :]
        
        laplacian = d2dx2 + d2dy2
        return laplacian
    
    def _find_energy_minima(self, energy_surface: np.ndarray, 
                           pc1_edges: np.ndarray, pc2_edges: np.ndarray,
                           min_separation: int = 3) -> List[Tuple[float, float, float]]:
        """
        Find local energy minima in the landscape
        
        Parameters:
        -----------
        energy_surface : np.ndarray
            2D energy landscape
        pc1_edges : np.ndarray
            PC1 bin edges
        pc2_edges : np.ndarray
            PC2 bin edges
        min_separation : int
            Minimum separation between minima (in bins)
            
        Returns:
        --------
        List[Tuple] : List of minima (pc1, pc2, energy)
        """
        from scipy.ndimage import minimum_filter
        
        # Find local minima using minimum filter
        local_min_mask = (energy_surface == minimum_filter(energy_surface, size=min_separation))
        
        # Get coordinates of minima
        min_coords = np.where(local_min_mask)
        
        minima = []
        for i, j in zip(min_coords[0], min_coords[1]):
            # Convert bin indices to PC coordinates
            pc1_coord = (pc1_edges[j] + pc1_edges[j+1]) / 2 if j < len(pc1_edges)-1 else pc1_edges[j]
            pc2_coord = (pc2_edges[i] + pc2_edges[i+1]) / 2 if i < len(pc2_edges)-1 else pc2_edges[i]
            energy = energy_surface[i, j]
            
            minima.append((pc1_coord, pc2_coord, energy))
        
        # Sort by energy (lowest first)
        minima.sort(key=lambda x: x[2])
        
        # Filter out minima that are too close to each other
        filtered_minima = []
        for minimum in minima:
            pc1_min, pc2_min, energy_min = minimum
            
            # Check if this minimum is far enough from existing ones
            too_close = False
            for existing_min in filtered_minima:
                pc1_exist, pc2_exist, _ = existing_min
                distance = np.sqrt((pc1_min - pc1_exist)**2 + (pc2_min - pc2_exist)**2)
                
                # Convert to bin units for comparison
                pc1_bin_size = (pc1_edges[-1] - pc1_edges[0]) / len(pc1_edges)
                pc2_bin_size = (pc2_edges[-1] - pc2_edges[0]) / len(pc2_edges)
                distance_bins = distance / np.sqrt(pc1_bin_size**2 + pc2_bin_size**2)
                
                if distance_bins < min_separation:
                    too_close = True
                    break
            
            if not too_close:
                filtered_minima.append(minimum)
        
        return filtered_minima[:10]  # Return top 10 minima
    
    def _analyze_energy_barriers(self, energy_surface: np.ndarray, 
                               minima: List[Tuple[float, float, float]],
                               pc1_edges: np.ndarray, pc2_edges: np.ndarray) -> List[Dict[str, Any]]:
        """
        Analyze energy barriers between minima
        
        Parameters:
        -----------
        energy_surface : np.ndarray
            2D energy landscape
        minima : List[Tuple]
            List of energy minima
        pc1_edges : np.ndarray
            PC1 bin edges  
        pc2_edges : np.ndarray
            PC2 bin edges
            
        Returns:
        --------
        List[Dict] : Energy barrier information
        """
        if len(minima) < 2:
            return []
        
        barriers = []
        
        # Analyze barriers between the lowest energy minima
        for i in range(min(len(minima), 5)):  # Top 5 minima
            for j in range(i+1, min(len(minima), 5)):
                min1 = minima[i]
                min2 = minima[j]
                
                # Find approximate barrier height
                pc1_path = np.linspace(min1[0], min2[0], 100)
                pc2_path = np.linspace(min1[1], min2[1], 100)
                
                # Convert to bin indices
                pc1_bin_size = (pc1_edges[-1] - pc1_edges[0]) / (len(pc1_edges) - 1)
                pc2_bin_size = (pc2_edges[-1] - pc2_edges[0]) / (len(pc2_edges) - 1)
                
                path_energies = []
                for pc1_val, pc2_val in zip(pc1_path, pc2_path):
                    # Find closest bin
                    i_bin = int((pc2_val - pc2_edges[0]) / pc2_bin_size)
                    j_bin = int((pc1_val - pc1_edges[0]) / pc1_bin_size)
                    
                    # Ensure within bounds
                    i_bin = max(0, min(i_bin, energy_surface.shape[0]-1))
                    j_bin = max(0, min(j_bin, energy_surface.shape[1]-1))
                    
                    path_energies.append(energy_surface[i_bin, j_bin])
                
                if path_energies:
                    barrier_height = max(path_energies) - min(min1[2], min2[2])
                    
                    barriers.append({
                        'minimum_1': min1,
                        'minimum_2': min2,
                        'barrier_height': barrier_height,
                        'min_energy_diff': abs(min1[2] - min2[2])
                    })
        
        # Sort by barrier height
        barriers.sort(key=lambda x: x['barrier_height'])
        
        return barriers
    
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
        from utils import timeout_handler
        
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
    
    def _advanced_community_detection(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """
        Advanced community detection using multiple algorithms
        
        Parameters:
        -----------
        G : nx.Graph
            Network graph
        metrics : NetworkMetrics
            Network metrics object to update
            
        Returns:
        --------
        NetworkMetrics : Updated metrics with community information
        """
        try:
            print(f"Computing communities using {self.config.community_method} method...")
            
            # Convert to numpy array for efficient computation
            nodes = list(G.nodes())
            node_to_idx = {node: i for i, node in enumerate(nodes)}
            
            # Create adjacency matrix
            adj_matrix = nx.adjacency_matrix(G, nodelist=nodes, weight='weight').toarray()
            
            if self.config.community_method == "leiden":
                communities, modularity = self._leiden_communities(G)
            elif self.config.community_method == "louvain":
                communities, modularity = self._louvain_communities(G)
            elif self.config.community_method == "spectral":
                communities, modularity = self._spectral_communities(G, adj_matrix, nodes)
            elif self.config.community_method == "hierarchical":
                communities, modularity = self._hierarchical_communities(G, adj_matrix, nodes)
            else:
                # Fallback to greedy modularity (NetworkX default)
                communities, modularity = self._greedy_modularity_communities(G)
            
            # Store detailed community information
            community_details = {
                'method': self.config.community_method,
                'n_communities': len(communities),
                'modularity': modularity,
                'communities': communities,
                'community_sizes': [len(comm) for comm in communities],
                'largest_community_size': max(len(comm) for comm in communities) if communities else 0,
                'smallest_community_size': min(len(comm) for comm in communities) if communities else 0
            }
            
            # Create node-to-community mapping
            node_assignments = {}
            for comm_id, community in enumerate(communities):
                for node in community:
                    node_assignments[node] = comm_id
            
            # Compute community statistics
            community_details['inter_community_edges'] = self._count_inter_community_edges(G, communities)
            community_details['intra_community_density'] = self._compute_intra_community_density(G, communities)
            community_details['inter_community_density'] = self._compute_inter_community_density(G, communities)
            
            # Store results
            metrics.communities = communities
            metrics.modularity = modularity
            metrics.communities_detailed = community_details
            metrics.community_node_assignments = node_assignments
            
            print(f"Found {len(communities)} communities with modularity = {modularity:.4f}")
            
        except Exception as e:
            print(f"Community detection failed: {e}")
            metrics.communities = []
            metrics.modularity = 0.0
            metrics.communities_detailed = {'method': 'failed', 'error': str(e)}
            metrics.community_node_assignments = {}
        
        return metrics
    
    def _leiden_communities(self, G: nx.Graph) -> Tuple[List[List[str]], float]:
        """Leiden community detection (requires python-igraph and leidenalg)"""
        try:
            import igraph as ig
            import leidenalg
            
            # Convert NetworkX to igraph
            edge_list = [(u, v, d['weight']) for u, v, d in G.edges(data=True)]
            vertices = list(G.nodes())
            vertex_map = {v: i for i, v in enumerate(vertices)}
            
            ig_graph = ig.Graph()
            ig_graph.add_vertices(len(vertices))
            ig_graph.add_edges([(vertex_map[u], vertex_map[v]) for u, v, _ in edge_list])
            ig_graph.es['weight'] = [w for _, _, w in edge_list]
            
            # Run Leiden algorithm
            partition = leidenalg.find_partition(ig_graph, leidenalg.ModularityVertexPartition, weights='weight')
            
            # Convert back to NetworkX format
            communities = []
            for community in partition:
                communities.append([vertices[i] for i in community])
            
            modularity = partition.modularity
            
            return communities, modularity
            
        except ImportError:
            print("Leiden algorithm requires igraph and leidenalg packages. Using Louvain instead.")
            return self._louvain_communities(G)
    
    def _louvain_communities(self, G: nx.Graph) -> Tuple[List[List[str]], float]:
        """Louvain community detection"""
        try:
            import networkx.algorithms.community as nx_comm
            
            # Use NetworkX's Louvain implementation
            communities_gen = nx_comm.louvain_communities(G, weight='weight', resolution=1.0)
            communities = [list(community) for community in communities_gen]
            
            # Compute modularity
            modularity = nx_comm.modularity(G, communities, weight='weight')
            
            return communities, modularity
            
        except Exception as e:
            print(f"Louvain failed: {e}, using greedy modularity")
            return self._greedy_modularity_communities(G)
    
    def _spectral_communities(self, G: nx.Graph, adj_matrix: np.ndarray, nodes: List[str]) -> Tuple[List[List[str]], float]:
        """Spectral clustering for community detection"""
        try:
            # Determine number of clusters using eigengap heuristic
            n_nodes = len(nodes)
            max_clusters = min(10, n_nodes // 5)  # Reasonable upper bound
            
            if max_clusters < 2:
                return [nodes], 0.0
            
            # Try different numbers of clusters and pick best modularity
            best_communities = None
            best_modularity = -1
            
            for n_clusters in range(2, max_clusters + 1):
                try:
                    # Spectral clustering
                    clustering = SpectralClustering(
                        n_clusters=n_clusters,
                        affinity='precomputed',
                        random_state=42
                    )
                    labels = clustering.fit_predict(adj_matrix)
                    
                    # Convert to community format
                    communities = defaultdict(list)
                    for i, label in enumerate(labels):
                        communities[label].append(nodes[i])
                    
                    community_list = list(communities.values())
                    
                    # Compute modularity
                    modularity = nx.community.modularity(G, community_list, weight='weight')
                    
                    if modularity > best_modularity:
                        best_modularity = modularity
                        best_communities = community_list
                        
                except Exception:
                    continue
            
            return best_communities or [nodes], best_modularity
            
        except Exception as e:
            print(f"Spectral clustering failed: {e}")
            return [nodes], 0.0
    
    def _hierarchical_communities(self, G: nx.Graph, adj_matrix: np.ndarray, nodes: List[str]) -> Tuple[List[List[str]], float]:
        """Hierarchical clustering for community detection"""
        try:
            # Convert adjacency matrix to distance matrix
            distance_matrix = 1.0 / (adj_matrix + 1e-10)  # Inverse of weights as distances
            np.fill_diagonal(distance_matrix, 0)
            
            # Hierarchical clustering
            condensed_distances = squareform(distance_matrix)
            linkage_matrix = linkage(condensed_distances, method='ward')
            
            # Try different numbers of clusters
            best_communities = None
            best_modularity = -1
            
            max_clusters = min(10, len(nodes) // 5)
            for n_clusters in range(2, max_clusters + 1):
                try:
                    labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
                    
                    # Convert to community format
                    communities = defaultdict(list)
                    for i, label in enumerate(labels):
                        communities[label].append(nodes[i])
                    
                    community_list = list(communities.values())
                    
                    # Compute modularity
                    modularity = nx.community.modularity(G, community_list, weight='weight')
                    
                    if modularity > best_modularity:
                        best_modularity = modularity
                        best_communities = community_list
                        
                except Exception:
                    continue
            
            return best_communities or [nodes], best_modularity
            
        except Exception as e:
            print(f"Hierarchical clustering failed: {e}")
            return [nodes], 0.0
    
    def _greedy_modularity_communities(self, G: nx.Graph) -> Tuple[List[List[str]], float]:
        """Fallback greedy modularity optimization"""
        try:
            communities_gen = nx.community.greedy_modularity_communities(G, weight='weight')
            communities = [list(community) for community in communities_gen]
            modularity = nx.community.modularity(G, communities, weight='weight')
            return communities, modularity
        except Exception as e:
            print(f"Greedy modularity failed: {e}")
            return [list(G.nodes())], 0.0
    
    def _count_inter_community_edges(self, G: nx.Graph, communities: List[List[str]]) -> int:
        """Count edges between different communities"""
        # Create community assignment mapping
        node_to_community = {}
        for i, community in enumerate(communities):
            for node in community:
                node_to_community[node] = i
        
        inter_edges = 0
        for u, v in G.edges():
            if node_to_community.get(u, -1) != node_to_community.get(v, -1):
                inter_edges += 1
        
        return inter_edges
    
    def _compute_intra_community_density(self, G: nx.Graph, communities: List[List[str]]) -> List[float]:
        """Compute density within each community"""
        densities = []
        for community in communities:
            if len(community) < 2:
                densities.append(0.0)
                continue
            
            subgraph = G.subgraph(community)
            density = nx.density(subgraph)
            densities.append(density)
        
        return densities
    
    def _compute_inter_community_density(self, G: nx.Graph, communities: List[List[str]]) -> float:
        """Compute density of edges between communities"""
        # Count total possible inter-community edges
        total_inter_possible = 0
        total_inter_actual = 0
        
        for i in range(len(communities)):
            for j in range(i + 1, len(communities)):
                comm1, comm2 = communities[i], communities[j]
                possible_edges = len(comm1) * len(comm2)
                total_inter_possible += possible_edges
                
                # Count actual edges
                actual_edges = 0
                for u in comm1:
                    for v in comm2:
                        if G.has_edge(u, v):
                            actual_edges += 1
                
                total_inter_actual += actual_edges
        
        return total_inter_actual / total_inter_possible if total_inter_possible > 0 else 0.0
    
    def _compute_path_metrics(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """
        Compute detailed path metrics for network analysis
        
        Parameters:
        -----------
        G : nx.Graph
            Network graph
        metrics : NetworkMetrics
            Network metrics object to update
            
        Returns:
        --------
        NetworkMetrics : Updated metrics with path information
        """
        try:
            print("Computing detailed path metrics...")
            
            path_metrics = {
                'network_efficiency': 0.0,
                'global_efficiency': 0.0,
                'local_efficiency': 0.0,
                'characteristic_path_length': 0.0,
                'network_diameter': 0.0,
                'node_eccentricities': {},
                'betweenness_distribution': {},
                'path_length_distribution': {},
                'shortest_paths_sample': {},
                'critical_paths': []
            }
            
            nodes = list(G.nodes())
            n_nodes = len(nodes)
            
            if n_nodes < 2:
                metrics.path_metrics = path_metrics
                return metrics
            
            # Compute all-pairs shortest paths (sample for large networks)
            if n_nodes > 500:
                # Sample nodes for efficiency
                sample_nodes = np.random.choice(nodes, min(100, n_nodes), replace=False)
                print(f"Sampling {len(sample_nodes)} nodes for path analysis")
            else:
                sample_nodes = nodes
            
            # Compute shortest path lengths
            path_lengths = []
            path_dict = {}
            
            for source in sample_nodes:
                try:
                    lengths = nx.single_source_shortest_path_length(G, source, weight='weight')
                    path_dict[source] = lengths
                    
                    for target, length in lengths.items():
                        if source != target and length > 0:
                            path_lengths.append(length)
                
                except nx.NetworkXNoPath:
                    continue
            
            if not path_lengths:
                metrics.path_metrics = path_metrics
                return metrics
            
            # Global path metrics
            path_metrics['characteristic_path_length'] = np.mean(path_lengths)
            path_metrics['network_diameter'] = max(path_lengths)
            path_metrics['global_efficiency'] = np.mean([1/length for length in path_lengths if length > 0])
            
            # Path length distribution
            unique_lengths, counts = np.unique(path_lengths, return_counts=True)
            path_metrics['path_length_distribution'] = {
                'lengths': unique_lengths.tolist(),
                'counts': counts.tolist(),
                'mean': np.mean(path_lengths),
                'std': np.std(path_lengths),
                'median': np.median(path_lengths)
            }
            
            # Node eccentricities (max distance from each node)
            for node in sample_nodes:
                if node in path_dict:
                    distances = [d for d in path_dict[node].values() if d > 0]
                    path_metrics['node_eccentricities'][node] = max(distances) if distances else 0
            
            # Local efficiency (average efficiency of node neighborhoods)
            local_efficiencies = []
            for node in sample_nodes:
                neighbors = list(G.neighbors(node))
                if len(neighbors) > 1:
                    subgraph = G.subgraph(neighbors)
                    if subgraph.number_of_edges() > 0:
                        try:
                            # Compute efficiency within neighborhood
                            neighbor_paths = []
                            for u in neighbors:
                                for v in neighbors:
                                    if u != v:
                                        try:
                                            length = nx.shortest_path_length(subgraph, u, v, weight='weight')
                                            neighbor_paths.append(1/length if length > 0 else 0)
                                        except nx.NetworkXNoPath:
                                            neighbor_paths.append(0)
                            
                            local_eff = np.mean(neighbor_paths) if neighbor_paths else 0
                            local_efficiencies.append(local_eff)
                        except:
                            local_efficiencies.append(0)
            
            path_metrics['local_efficiency'] = np.mean(local_efficiencies) if local_efficiencies else 0.0
            
            # Network efficiency
            path_metrics['network_efficiency'] = 1 / path_metrics['characteristic_path_length'] if path_metrics['characteristic_path_length'] > 0 else 0
            
            # Identify critical paths (paths through high-betweenness nodes)
            if metrics.betweenness_centrality:
                high_betweenness_nodes = sorted(metrics.betweenness_centrality.items(), 
                                              key=lambda x: x[1], reverse=True)[:5]
                
                critical_paths = []
                for node, betweenness in high_betweenness_nodes:
                    if node in sample_nodes and node in path_dict:
                        # Find paths that go through this high-betweenness node
                        node_paths = []
                        for target in path_dict[node]:
                            if path_dict[node][target] > 0:
                                try:
                                    path = nx.shortest_path(G, node, target, weight='weight')
                                    node_paths.append({
                                        'source': node,
                                        'target': target,
                                        'length': path_dict[node][target],
                                        'path': path,
                                        'hub_betweenness': betweenness
                                    })
                                except:
                                    continue
                        
                        # Keep top paths through this hub
                        node_paths.sort(key=lambda x: x['length'], reverse=True)
                        critical_paths.extend(node_paths[:3])  # Top 3 longest paths
                
                path_metrics['critical_paths'] = critical_paths[:10]  # Top 10 overall
            
            # Store sample paths for specific node pairs if specified
            if self.config.path_analysis_nodes:
                sample_paths = {}
                analysis_nodes = [n for n in self.config.path_analysis_nodes if n in G.nodes()]
                
                for i, source in enumerate(analysis_nodes):
                    for target in analysis_nodes[i+1:]:
                        try:
                            path = nx.shortest_path(G, source, target, weight='weight')
                            length = nx.shortest_path_length(G, source, target, weight='weight')
                            
                            sample_paths[f"{source}_to_{target}"] = {
                                'path': path,
                                'length': length,
                                'n_steps': len(path) - 1
                            }
                        except nx.NetworkXNoPath:
                            sample_paths[f"{source}_to_{target}"] = {
                                'path': None,
                                'length': float('inf'),
                                'n_steps': -1
                            }
                
                path_metrics['shortest_paths_sample'] = sample_paths
            
            metrics.path_metrics = path_metrics
            print(f"Path analysis complete: diameter = {path_metrics['network_diameter']:.2f}, "
                  f"avg_path = {path_metrics['characteristic_path_length']:.2f}")
            
        except Exception as e:
            print(f"Path metrics computation failed: {e}")
            metrics.path_metrics = {'error': str(e)}
        
        return metrics
    
    def _analyze_allosteric_pathways(self, G: nx.Graph, metrics: NetworkMetrics, 
                                   simulation: MDSimulation) -> NetworkMetrics:
        """
        Analyze allosteric communication pathways in the protein network
        
        Parameters:
        -----------
        G : nx.Graph
            Network graph
        metrics : NetworkMetrics
            Network metrics object to update
        simulation : MDSimulation
            Simulation object with residue information
            
        Returns:
        --------
        NetworkMetrics : Updated metrics with allosteric pathway information
        """
        try:
            print("Analyzing allosteric communication pathways...")
            
            allosteric_results = {
                'source_nodes': [],
                'target_nodes': [], 
                'pathways': [],
                'communication_efficiency': {},
                'pathway_redundancy': {},
                'allosteric_hotspots': [],
                'functional_regions': {},
                'pathway_robustness': {}
            }
            
            # Define source and target nodes for allosteric analysis
            source_nodes = self._identify_allosteric_sources(G, metrics, simulation)
            target_nodes = self._identify_allosteric_targets(G, metrics, simulation)
            
            # Use user-defined nodes if provided
            if self.config.allosteric_source_nodes:
                user_sources = [n for n in self.config.allosteric_source_nodes if n in G.nodes()]
                if user_sources:
                    source_nodes = user_sources
            
            if self.config.allosteric_target_nodes:
                user_targets = [n for n in self.config.allosteric_target_nodes if n in G.nodes()]
                if user_targets:
                    target_nodes = user_targets
            
            allosteric_results['source_nodes'] = source_nodes
            allosteric_results['target_nodes'] = target_nodes
            
            if not source_nodes or not target_nodes:
                print("Warning: No source or target nodes identified for allosteric analysis")
                metrics.allosteric_pathways = allosteric_results
                return metrics
            
            # Compute allosteric pathways
            pathways = []
            communication_scores = {}
            
            for source in source_nodes:
                for target in target_nodes:
                    if source != target:
                        pathway_info = self._compute_allosteric_pathway(G, source, target, metrics)
                        if pathway_info:
                            pathways.append(pathway_info)
                            
                            # Store communication efficiency
                            pair_key = f"{source}_to_{target}"
                            communication_scores[pair_key] = pathway_info.get('efficiency', 0.0)
            
            allosteric_results['pathways'] = pathways
            allosteric_results['communication_efficiency'] = communication_scores
            
            # Identify allosteric hotspots (nodes frequently on pathways)
            hotspots = self._identify_allosteric_hotspots(pathways, G)
            allosteric_results['allosteric_hotspots'] = hotspots
            
            # Analyze pathway redundancy
            redundancy = self._analyze_pathway_redundancy(pathways, G)
            allosteric_results['pathway_redundancy'] = redundancy
            
            # Identify functional regions based on communities and pathways
            functional_regions = self._identify_functional_regions(pathways, metrics.communities or [])
            allosteric_results['functional_regions'] = functional_regions
            
            # Assess pathway robustness to node removal
            robustness = self._assess_pathway_robustness(G, pathways[:5])  # Top 5 pathways
            allosteric_results['pathway_robustness'] = robustness
            
            metrics.allosteric_pathways = allosteric_results
            
            print(f"Allosteric analysis complete: {len(pathways)} pathways identified, "
                  f"{len(hotspots)} hotspot residues")
            
        except Exception as e:
            print(f"Allosteric pathway analysis failed: {e}")
            metrics.allosteric_pathways = {'error': str(e)}
        
        return metrics
    
    def _identify_allosteric_sources(self, G: nx.Graph, metrics: NetworkMetrics, 
                                   simulation: MDSimulation) -> List[str]:
        """Identify potential allosteric source nodes (e.g., binding sites, flaps)"""
        sources = []
        
        # Use high-degree centrality nodes as potential sources
        if metrics.degree_centrality:
            high_degree = sorted(metrics.degree_centrality.items(), key=lambda x: x[1], reverse=True)
            sources.extend([node for node, _ in high_degree[:10]])
        
        # For HIV protease, add known functional regions
        if hasattr(simulation, 'unique_residue_keys'):
            # Flap regions (residues 45-55 in HIV protease)
            flap_residues = [res for res in simulation.unique_residue_keys 
                           if any(f'_{i}' in res for i in range(45, 56))]
            sources.extend(flap_residues)
            
            # Substrate binding region
            binding_residues = [res for res in simulation.unique_residue_keys 
                              if any(f'_{i}' in res for i in [8, 23, 25, 27, 29, 30])]
            sources.extend(binding_residues)
        
        return list(set(sources))[:15]  # Limit to top 15 sources
    
    def _identify_allosteric_targets(self, G: nx.Graph, metrics: NetworkMetrics, 
                                   simulation: MDSimulation) -> List[str]:
        """Identify potential allosteric target nodes (e.g., active site)"""
        targets = []
        
        # Use high-betweenness centrality nodes as potential targets
        if metrics.betweenness_centrality:
            high_betweenness = sorted(metrics.betweenness_centrality.items(), key=lambda x: x[1], reverse=True)
            targets.extend([node for node, _ in high_betweenness[:10]])
        
        # For HIV protease, add known functional regions
        if hasattr(simulation, 'unique_residue_keys'):
            # Active site (catalytic residues)
            active_site = [res for res in simulation.unique_residue_keys 
                          if any(f'_{i}' in res for i in [25, 27])]
            targets.extend(active_site)
            
            # Interface residues
            interface_residues = [res for res in simulation.unique_residue_keys 
                                if any(f'_{i}' in res for i in [80, 82, 84, 90])]
            targets.extend(interface_residues)
        
        return list(set(targets))[:15]  # Limit to top 15 targets
    
    def _compute_allosteric_pathway(self, G: nx.Graph, source: str, target: str, 
                                  metrics: NetworkMetrics) -> Dict[str, Any]:
        """Compute detailed allosteric pathway between source and target"""
        try:
            # Find shortest path
            shortest_path = nx.shortest_path(G, source, target, weight='weight')
            shortest_length = nx.shortest_path_length(G, source, target, weight='weight')
            
            # Analyze pathway characteristics
            pathway_betweenness = []
            if metrics.betweenness_centrality:
                pathway_betweenness = [metrics.betweenness_centrality.get(node, 0) for node in shortest_path]
            
            # Communication efficiency (inverse of path length)
            efficiency = 1.0 / shortest_length if shortest_length > 0 else 0.0
            
            pathway_info = {
                'source': source,
                'target': target,
                'shortest_path': shortest_path,
                'shortest_length': shortest_length,
                'path_nodes': len(shortest_path),
                'efficiency': efficiency,
                'pathway_betweenness': pathway_betweenness,
                'avg_pathway_betweenness': np.mean(pathway_betweenness) if pathway_betweenness else 0.0,
                'bottleneck_nodes': [node for node, bet in zip(shortest_path, pathway_betweenness) 
                                   if bet > np.mean(pathway_betweenness)] if pathway_betweenness else []
            }
            
            return pathway_info
            
        except nx.NetworkXNoPath:
            return None
        except Exception as e:
            print(f"Pathway computation failed for {source} -> {target}: {e}")
            return None
    
    def _identify_allosteric_hotspots(self, pathways: List[Dict[str, Any]], G: nx.Graph) -> List[Dict[str, Any]]:
        """Identify nodes that frequently appear in allosteric pathways"""
        node_frequency = defaultdict(int)
        node_pathway_info = defaultdict(list)
        
        # Count how often each node appears in pathways
        for pathway in pathways:
            if 'shortest_path' in pathway:
                for node in pathway['shortest_path']:
                    node_frequency[node] += 1
                    node_pathway_info[node].append({
                        'source': pathway['source'],
                        'target': pathway['target'],
                        'efficiency': pathway.get('efficiency', 0.0)
                    })
        
        # Create hotspot information
        hotspots = []
        for node, frequency in node_frequency.items():
            if frequency >= 2:  # Appears in at least 2 pathways
                avg_efficiency = np.mean([info['efficiency'] for info in node_pathway_info[node]])
                
                hotspots.append({
                    'node': node,
                    'frequency': frequency,
                    'pathway_count': len(node_pathway_info[node]),
                    'avg_communication_efficiency': avg_efficiency,
                    'hotspot_score': frequency * avg_efficiency
                })
        
        # Sort by hotspot score
        hotspots.sort(key=lambda x: x['hotspot_score'], reverse=True)
        
        return hotspots[:20]  # Top 20 hotspots
    
    def _analyze_pathway_redundancy(self, pathways: List[Dict[str, Any]], G: nx.Graph) -> Dict[str, Any]:
        """Analyze redundancy in allosteric communication"""
        redundancy_info = {
            'total_pathways': len(pathways),
            'unique_source_target_pairs': 0,
            'communication_robustness': 0.0
        }
        
        # Group pathways by source-target pairs
        pair_pathways = defaultdict(list)
        for pathway in pathways:
            pair_key = f"{pathway['source']}_to_{pathway['target']}"
            pair_pathways[pair_key].append(pathway)
        
        redundancy_info['unique_source_target_pairs'] = len(pair_pathways)
        redundancy_info['communication_robustness'] = len(pair_pathways) / len(pathways) if pathways else 0.0
        
        return redundancy_info
    
    def _identify_functional_regions(self, pathways: List[Dict[str, Any]], 
                                   communities: List[List[str]]) -> Dict[str, Any]:
        """Identify functional regions based on pathway and community analysis"""
        functional_regions = {
            'pathway_communities': [],
            'inter_community_bridges': []
        }
        
        if not communities or not pathways:
            return functional_regions
        
        # Create node-to-community mapping
        node_to_comm = {}
        for i, community in enumerate(communities):
            for node in community:
                node_to_comm[node] = i
        
        # Find nodes that bridge communities in pathways
        bridge_nodes = set()
        for pathway in pathways:
            if 'shortest_path' in pathway:
                for i in range(len(pathway['shortest_path']) - 1):
                    node = pathway['shortest_path'][i]
                    next_node = pathway['shortest_path'][i + 1]
                    
                    node_comm = node_to_comm.get(node)
                    next_comm = node_to_comm.get(next_node)
                    
                    if node_comm is not None and next_comm is not None and node_comm != next_comm:
                        bridge_nodes.add(node)
                        bridge_nodes.add(next_node)
        
        functional_regions['inter_community_bridges'] = list(bridge_nodes)
        
        return functional_regions
    
    def _assess_pathway_robustness(self, G: nx.Graph, pathways: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Assess robustness of pathways to node removal"""
        robustness_info = {
            'critical_nodes': [],
            'network_resilience': 0.0
        }
        
        # Simplified robustness analysis
        critical_counts = defaultdict(int)
        
        for pathway in pathways:
            if 'shortest_path' not in pathway or len(pathway['shortest_path']) < 3:
                continue
            
            # Mark intermediate nodes as critical
            for node in pathway['shortest_path'][1:-1]:  # Skip source and target
                critical_counts[node] += 1
        
        robustness_info['critical_nodes'] = sorted(critical_counts.items(), 
                                                 key=lambda x: x[1], reverse=True)[:10]
        
        # Simple resilience measure
        total_nodes = len(G.nodes())
        critical_node_count = len([node for node, count in critical_counts.items() if count >= 2])
        robustness_info['network_resilience'] = 1.0 - (critical_node_count / total_nodes) if total_nodes > 0 else 0.0
        
        return robustness_info
    
    def _analyze_network_robustness(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """
        Analyze network robustness to node/edge removal
        
        Parameters:
        -----------
        G : nx.Graph
            Network graph
        metrics : NetworkMetrics
            Network metrics object to update
            
        Returns:
        --------
        NetworkMetrics : Updated metrics with robustness information
        """
        try:
            print("Analyzing network robustness...")
            
            robustness_metrics = {
                'node_attack_robustness': 0.0,
                'random_failure_robustness': 0.0,
                'edge_attack_robustness': 0.0,
                'giant_component_size': 0.0,
                'fragmentation_point': None,
                'critical_nodes': [],
                'critical_edges': []
            }
            
            # Giant component size
            if nx.is_connected(G):
                robustness_metrics['giant_component_size'] = 1.0
            else:
                components = list(nx.connected_components(G))
                largest_component_size = max(len(comp) for comp in components) if components else 0
                robustness_metrics['giant_component_size'] = largest_component_size / len(G.nodes()) if G.nodes() else 0.0
            
            # Node attack robustness (remove highest centrality nodes)
            if metrics.betweenness_centrality and len(G.nodes()) > 5:
                G_test = G.copy()
                high_centrality_nodes = sorted(metrics.betweenness_centrality.items(), 
                                             key=lambda x: x[1], reverse=True)
                
                connectivity_loss = []
                nodes_removed = 0
                max_removals = min(10, len(G.nodes()) // 4)  # Remove up to 25% or 10 nodes
                
                for node, _ in high_centrality_nodes[:max_removals]:
                    G_test.remove_node(node)
                    nodes_removed += 1
                    
                    if nx.is_connected(G_test):
                        connectivity_loss.append(0)  # Still connected
                    else:
                        components = list(nx.connected_components(G_test))
                        largest_comp = max(len(comp) for comp in components) if components else 0
                        loss = 1.0 - (largest_comp / len(G.nodes()))
                        connectivity_loss.append(loss)
                        
                        if loss > 0.5 and robustness_metrics['fragmentation_point'] is None:
                            robustness_metrics['fragmentation_point'] = nodes_removed
                
                robustness_metrics['node_attack_robustness'] = 1.0 - np.mean(connectivity_loss) if connectivity_loss else 1.0
                robustness_metrics['critical_nodes'] = [node for node, _ in high_centrality_nodes[:5]]
            
            # Random failure robustness (simplified)
            if len(G.nodes()) > 10:
                nodes = list(G.nodes())
                random_nodes = np.random.choice(nodes, min(5, len(nodes) // 5), replace=False)
                
                G_random = G.copy()
                for node in random_nodes:
                    G_random.remove_node(node)
                
                if nx.is_connected(G_random):
                    robustness_metrics['random_failure_robustness'] = 1.0
                else:
                    components = list(nx.connected_components(G_random))
                    largest_comp = max(len(comp) for comp in components) if components else 0
                    robustness_metrics['random_failure_robustness'] = largest_comp / len(G.nodes()) if G.nodes() else 0.0
            
            # Edge attack robustness (remove high-betweenness edges)
            edge_betweenness = nx.edge_betweenness_centrality(G, weight='weight', k=min(50, len(G.nodes())))
            if edge_betweenness:
                high_bet_edges = sorted(edge_betweenness.items(), key=lambda x: x[1], reverse=True)[:5]
                robustness_metrics['critical_edges'] = [edge for edge, _ in high_bet_edges]
                
                G_edge_test = G.copy()
                for edge, _ in high_bet_edges[:3]:  # Remove top 3 edges
                    if G_edge_test.has_edge(*edge):
                        G_edge_test.remove_edge(*edge)
                
                if nx.is_connected(G_edge_test):
                    robustness_metrics['edge_attack_robustness'] = 1.0
                else:
                    components = list(nx.connected_components(G_edge_test))
                    largest_comp = max(len(comp) for comp in components) if components else 0
                    robustness_metrics['edge_attack_robustness'] = largest_comp / len(G.nodes()) if G.nodes() else 0.0
            
            metrics.network_robustness = robustness_metrics
            print(f"Network robustness analysis complete: "
                  f"node_attack={robustness_metrics['node_attack_robustness']:.3f}, "
                  f"giant_component={robustness_metrics['giant_component_size']:.3f}")
            
        except Exception as e:
            print(f"Network robustness analysis failed: {e}")
            metrics.network_robustness = {'error': str(e)}
        
        return metrics
    
    def _compute_centrality_significance(self, G: nx.Graph, metrics: NetworkMetrics) -> NetworkMetrics:
        """
        Compute statistical significance of centrality measures using z-scores
        
        Parameters:
        -----------
        G : nx.Graph
            Network graph
        metrics : NetworkMetrics
            Network metrics object to update
            
        Returns:
        --------
        NetworkMetrics : Updated metrics with centrality significance scores
        """
        try:
            print("Computing centrality significance scores...")
            
            z_scores = {
                'betweenness_z': {},
                'closeness_z': {},
                'degree_z': {},
                'eigenvector_z': {}
            }
            
            # Calculate z-scores for each centrality measure
            centrality_measures = {
                'betweenness': metrics.betweenness_centrality or {},
                'closeness': metrics.closeness_centrality or {},
                'degree': metrics.degree_centrality or {},
                'eigenvector': metrics.eigenvector_centrality or {}
            }
            
            for measure_name, centrality_dict in centrality_measures.items():
                if centrality_dict:
                    values = list(centrality_dict.values())
                    
                    if len(values) > 1:
                        mean_val = np.mean(values)
                        std_val = np.std(values)
                        
                        if std_val > 0:
                            for node, value in centrality_dict.items():
                                z_score = (value - mean_val) / std_val
                                z_scores[f"{measure_name}_z"][node] = z_score
            
            # Identify statistically significant nodes (|z| > 2)
            significant_nodes = {
                'highly_significant': set(),  # |z| > 2.5
                'significant': set(),         # 1.96 < |z| <= 2.5
                'notable': set()             # 1.0 < |z| <= 1.96
            }
            
            for measure_name, z_dict in z_scores.items():
                for node, z_score in z_dict.items():
                    abs_z = abs(z_score)
                    
                    if abs_z > 2.5:
                        significant_nodes['highly_significant'].add(node)
                    elif abs_z > 1.96:
                        significant_nodes['significant'].add(node)
                    elif abs_z > 1.0:
                        significant_nodes['notable'].add(node)
            
            # Add significance categories to z_scores
            z_scores['significant_nodes'] = {
                'highly_significant': list(significant_nodes['highly_significant']),
                'significant': list(significant_nodes['significant']),
                'notable': list(significant_nodes['notable'])
            }
            
            # Summary statistics
            z_scores['summary'] = {
                'n_highly_significant': len(significant_nodes['highly_significant']),
                'n_significant': len(significant_nodes['significant']),
                'n_notable': len(significant_nodes['notable']),
                'total_significant': len(significant_nodes['highly_significant']) + len(significant_nodes['significant'])
            }
            
            metrics.centrality_z_scores = z_scores
            
            total_sig = z_scores['summary']['total_significant']
            print(f"Centrality significance analysis complete: {total_sig} statistically significant nodes")
            
        except Exception as e:
            print(f"Centrality significance analysis failed: {e}")
            metrics.centrality_z_scores = {'error': str(e)}
        
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
        
        # Save dynamic analysis results
        if hasattr(simulation, 'dynamic_analysis') and simulation.dynamic_analysis:
            self._save_dynamic_analysis_data(simulation, prefix)
        
        # Save advanced network analysis results
        self._save_advanced_network_data(simulation, prefix)
        
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
    
    def _save_dynamic_analysis_data(self, simulation: MDSimulation, prefix: str):
        """Save dynamic analysis results (DCCM and PCA)"""
        if not hasattr(simulation, 'dynamic_analysis') or not simulation.dynamic_analysis:
            return
            
        dynamic_data = simulation.dynamic_analysis
        
        # Save DCCM matrix
        if 'dccm_matrix' in dynamic_data and dynamic_data['dccm_matrix'] is not None:
            np.save(self.output_dir / f"{prefix}_dccm_matrix.npy", dynamic_data['dccm_matrix'])
            
            # Save as CSV if not too large
            if dynamic_data['dccm_matrix'].size < 1000000:
                np.savetxt(
                    self.output_dir / f"{prefix}_dccm_matrix.csv",
                    dynamic_data['dccm_matrix'],
                    delimiter=",",
                    fmt="%.6f",
                    comments=""
                )
            
            # Create DCCM visualization
            self._create_dccm_visualization(dynamic_data['dccm_matrix'], prefix)
        
        # Save PCA results
        if 'pca_eigenvalues' in dynamic_data and dynamic_data['pca_eigenvalues'] is not None:
            np.save(self.output_dir / f"{prefix}_pca_eigenvalues.npy", dynamic_data['pca_eigenvalues'])
            np.save(self.output_dir / f"{prefix}_pca_eigenvectors.npy", dynamic_data['pca_eigenvectors'])
            np.save(self.output_dir / f"{prefix}_pca_variance_explained.npy", dynamic_data['pca_variance_explained'])
            
            if 'principal_components' in dynamic_data and dynamic_data['principal_components'] is not None:
                np.save(self.output_dir / f"{prefix}_principal_components.npy", dynamic_data['principal_components'])
            
            # Save PCA summary
            with open(self.output_dir / f"{prefix}_pca_summary.txt", 'w') as f:
                f.write("Principal Component Analysis Summary\n")
                f.write("=" * 40 + "\n")
                f.write(f"Number of components: {dynamic_data.get('pca_n_components', 'N/A')}\n")
                f.write(f"Total variance: {dynamic_data.get('pca_total_variance', 'N/A'):.6f}\n")
                
                if 'pca_cumulative_variance' in dynamic_data:
                    f.write(f"Variance explained by first 5 PCs:\n")
                    for i, var_exp in enumerate(dynamic_data['pca_variance_explained'][:5]):
                        f.write(f"  PC{i+1}: {var_exp*100:.2f}%\n")
                    
                    f.write(f"\nCumulative variance explained:\n")
                    for i, cum_var in enumerate(dynamic_data['pca_cumulative_variance'][:5]):
                        f.write(f"  PC1-{i+1}: {cum_var*100:.2f}%\n")
            
            # Create PCA visualization
            self._create_pca_visualization(dynamic_data, prefix)
        
        # Save energy landscape results
        if 'energy_landscape' in dynamic_data and dynamic_data['energy_landscape'] is not None:
            # Save landscape data
            np.save(self.output_dir / f"{prefix}_energy_landscape.npy", dynamic_data['energy_landscape'])
            np.save(self.output_dir / f"{prefix}_landscape_pc1_bins.npy", dynamic_data['landscape_pc1_bins'])
            np.save(self.output_dir / f"{prefix}_landscape_pc2_bins.npy", dynamic_data['landscape_pc2_bins'])
            np.save(self.output_dir / f"{prefix}_landscape_gradient_x.npy", dynamic_data['landscape_gradient_x'])
            np.save(self.output_dir / f"{prefix}_landscape_gradient_y.npy", dynamic_data['landscape_gradient_y'])
            np.save(self.output_dir / f"{prefix}_landscape_laplacian.npy", dynamic_data['landscape_laplacian'])
            
            # Save landscape summary
            with open(self.output_dir / f"{prefix}_landscape_summary.txt", 'w') as f:
                f.write("Energy Landscape Analysis Summary\n")
                f.write("=" * 40 + "\n")
                f.write(f"Temperature: {dynamic_data.get('landscape_temperature', 'N/A')} K\n")
                f.write(f"Total frames analyzed: {dynamic_data.get('landscape_total_frames', 'N/A')}\n")
                
                energy_range = dynamic_data.get('landscape_energy_range', [0, 0])
                f.write(f"Energy range: {energy_range[0]:.1f} to {energy_range[1]:.1f} kJ/mol\n")
                
                # Minima information
                minima = dynamic_data.get('landscape_minima', [])
                f.write(f"\nNumber of energy minima found: {len(minima)}\n")
                if minima:
                    f.write("Top 5 energy minima (PC1, PC2, Energy):\n")
                    for i, (pc1, pc2, energy) in enumerate(minima[:5]):
                        f.write(f"  {i+1}: PC1={pc1:.3f}, PC2={pc2:.3f}, E={energy:.1f} kJ/mol\n")
                
                # Barrier information
                barriers = dynamic_data.get('landscape_barriers', [])
                f.write(f"\nNumber of energy barriers analyzed: {len(barriers)}\n")
                if barriers:
                    f.write("Top 3 energy barriers:\n")
                    for i, barrier in enumerate(barriers[:3]):
                        height = barrier['barrier_height']
                        f.write(f"  {i+1}: Barrier height = {height:.1f} kJ/mol\n")
            
            # Create landscape visualizations
            self._create_landscape_visualization(dynamic_data, prefix)
        
        # Save dynamic analysis summary
        with open(self.output_dir / f"{prefix}_dynamic_summary.txt", 'w') as f:
            f.write("Dynamic Analysis Summary\n")
            f.write("=" * 30 + "\n")
            
            if 'dccm_matrix' in dynamic_data:
                dccm = dynamic_data['dccm_matrix']
                f.write(f"DCCM computed: Yes\n")
                f.write(f"DCCM atoms: {dynamic_data.get('dccm_atoms', 'N/A')}\n")
                f.write(f"DCCM matrix size: {dccm.shape}\n")
                f.write(f"DCCM range: {np.min(dccm):.3f} to {np.max(dccm):.3f}\n")
                f.write(f"Average correlation: {np.mean(dccm):.3f}\n")
            else:
                f.write(f"DCCM computed: No\n")
            
            if 'pca_eigenvalues' in dynamic_data:
                f.write(f"PCA computed: Yes\n")
                f.write(f"PCA components: {dynamic_data.get('pca_n_components', 'N/A')}\n")
                if 'pca_cumulative_variance' in dynamic_data and len(dynamic_data['pca_cumulative_variance']) > 0:
                    f.write(f"First PC explains: {dynamic_data['pca_variance_explained'][0]*100:.1f}% variance\n")
                    f.write(f"First 3 PCs explain: {dynamic_data['pca_cumulative_variance'][2]*100:.1f}% variance\n")
            else:
                f.write(f"PCA computed: No\n")
                
            if 'energy_landscape' in dynamic_data:
                f.write(f"Energy landscape computed: Yes\n")
                energy_range = dynamic_data.get('landscape_energy_range', [0, 0])
                f.write(f"Energy range: {energy_range[1] - energy_range[0]:.1f} kJ/mol\n")
                minima_count = len(dynamic_data.get('landscape_minima', []))
                f.write(f"Energy minima found: {minima_count}\n")
            else:
                f.write(f"Energy landscape computed: No\n")
    
    def _save_advanced_network_data(self, simulation: MDSimulation, prefix: str):
        """Save advanced network analysis results"""
        if not simulation.network_metrics:
            return
            
        metrics = simulation.network_metrics
        
        # Save community detection results
        if hasattr(metrics, 'communities_detailed') and metrics.communities_detailed:
            with open(self.output_dir / f"{prefix}_community_analysis.json", 'w') as f:
                json.dump(metrics.communities_detailed, f, indent=2, default=str)
            
            # Save community assignments
            if hasattr(metrics, 'community_node_assignments') and metrics.community_node_assignments:
                with open(self.output_dir / f"{prefix}_community_assignments.csv", 'w') as f:
                    f.write("node,community_id\n")
                    for node, comm_id in metrics.community_node_assignments.items():
                        f.write(f"{node},{comm_id}\n")
            
            # Save community summary
            with open(self.output_dir / f"{prefix}_community_summary.txt", 'w') as f:
                f.write("Community Detection Analysis Summary\n")
                f.write("=" * 45 + "\n")
                f.write(f"Method: {metrics.communities_detailed.get('method', 'N/A')}\n")
                f.write(f"Number of communities: {metrics.communities_detailed.get('n_communities', 0)}\n")
                f.write(f"Modularity: {metrics.communities_detailed.get('modularity', 0.0):.4f}\n")
                
                sizes = metrics.communities_detailed.get('community_sizes', [])
                if sizes:
                    f.write(f"Largest community: {max(sizes)} nodes\n")
                    f.write(f"Smallest community: {min(sizes)} nodes\n")
                    f.write(f"Average community size: {np.mean(sizes):.1f} nodes\n")
                
                inter_edges = metrics.communities_detailed.get('inter_community_edges', 0)
                f.write(f"Inter-community edges: {inter_edges}\n")
        
        # Save path metrics
        if hasattr(metrics, 'path_metrics') and metrics.path_metrics:
            with open(self.output_dir / f"{prefix}_path_metrics.json", 'w') as f:
                json.dump(metrics.path_metrics, f, indent=2, default=str)
            
            with open(self.output_dir / f"{prefix}_path_summary.txt", 'w') as f:
                f.write("Path Analysis Summary\n")
                f.write("=" * 25 + "\n")
                pm = metrics.path_metrics
                f.write(f"Characteristic path length: {pm.get('characteristic_path_length', 0):.3f}\n")
                f.write(f"Network diameter: {pm.get('network_diameter', 0):.3f}\n")
                f.write(f"Global efficiency: {pm.get('global_efficiency', 0):.4f}\n")
                f.write(f"Local efficiency: {pm.get('local_efficiency', 0):.4f}\n")
                f.write(f"Network efficiency: {pm.get('network_efficiency', 0):.4f}\n")
                
                # Critical paths
                critical_paths = pm.get('critical_paths', [])
                if critical_paths:
                    f.write(f"\nTop {min(5, len(critical_paths))} Critical Paths:\n")
                    for i, path_info in enumerate(critical_paths[:5]):
                        f.write(f"  {i+1}. {path_info['source']} → {path_info['target']} "
                               f"(length: {path_info['length']:.2f})\n")
        
        # Save allosteric pathway analysis
        if hasattr(metrics, 'allosteric_pathways') and metrics.allosteric_pathways:
            with open(self.output_dir / f"{prefix}_allosteric_analysis.json", 'w') as f:
                json.dump(metrics.allosteric_pathways, f, indent=2, default=str)
            
            with open(self.output_dir / f"{prefix}_allosteric_summary.txt", 'w') as f:
                f.write("Allosteric Pathway Analysis Summary\n")
                f.write("=" * 40 + "\n")
                ap = metrics.allosteric_pathways
                
                sources = ap.get('source_nodes', [])
                targets = ap.get('target_nodes', [])
                pathways = ap.get('pathways', [])
                hotspots = ap.get('allosteric_hotspots', [])
                
                f.write(f"Source nodes identified: {len(sources)}\n")
                f.write(f"Target nodes identified: {len(targets)}\n")
                f.write(f"Allosteric pathways found: {len(pathways)}\n")
                f.write(f"Allosteric hotspots: {len(hotspots)}\n")
                
                if sources:
                    f.write(f"\nTop source nodes: {', '.join(sources[:5])}\n")
                if targets:
                    f.write(f"Top target nodes: {', '.join(targets[:5])}\n")
                
                if hotspots:
                    f.write(f"\nTop 5 Allosteric Hotspots:\n")
                    for i, hotspot in enumerate(hotspots[:5]):
                        f.write(f"  {i+1}. {hotspot['node']} "
                               f"(score: {hotspot['hotspot_score']:.3f}, "
                               f"frequency: {hotspot['frequency']})\n")
                
                # Communication efficiency
                comm_eff = ap.get('communication_efficiency', {})
                if comm_eff:
                    avg_efficiency = np.mean(list(comm_eff.values()))
                    f.write(f"\nAverage communication efficiency: {avg_efficiency:.4f}\n")
                    
                    # Top communication pairs
                    top_pairs = sorted(comm_eff.items(), key=lambda x: x[1], reverse=True)[:3]
                    f.write("Most efficient communication pairs:\n")
                    for pair, eff in top_pairs:
                        f.write(f"  {pair}: {eff:.4f}\n")
            
            # Save hotspots as CSV for easy analysis
            if hotspots:
                with open(self.output_dir / f"{prefix}_allosteric_hotspots.csv", 'w') as f:
                    f.write("node,hotspot_score,frequency,pathway_count,avg_efficiency\n")
                    for hotspot in hotspots:
                        f.write(f"{hotspot['node']},{hotspot['hotspot_score']:.4f},"
                               f"{hotspot['frequency']},{hotspot['pathway_count']},"
                               f"{hotspot['avg_communication_efficiency']:.4f}\n")
        
        # Save network robustness analysis
        if hasattr(metrics, 'network_robustness') and metrics.network_robustness:
            with open(self.output_dir / f"{prefix}_robustness_analysis.json", 'w') as f:
                json.dump(metrics.network_robustness, f, indent=2, default=str)
            
            with open(self.output_dir / f"{prefix}_robustness_summary.txt", 'w') as f:
                f.write("Network Robustness Analysis Summary\n")
                f.write("=" * 40 + "\n")
                nr = metrics.network_robustness
                
                f.write(f"Giant component size: {nr.get('giant_component_size', 0):.3f}\n")
                f.write(f"Node attack robustness: {nr.get('node_attack_robustness', 0):.3f}\n")
                f.write(f"Random failure robustness: {nr.get('random_failure_robustness', 0):.3f}\n")
                f.write(f"Edge attack robustness: {nr.get('edge_attack_robustness', 0):.3f}\n")
                
                frag_point = nr.get('fragmentation_point')
                if frag_point:
                    f.write(f"Network fragmentation point: {frag_point} nodes removed\n")
                
                critical_nodes = nr.get('critical_nodes', [])
                if critical_nodes:
                    f.write(f"\nTop 5 Critical Nodes:\n")
                    for i, node in enumerate(critical_nodes[:5]):
                        f.write(f"  {i+1}. {node}\n")
        
        # Save centrality significance scores
        if hasattr(metrics, 'centrality_z_scores') and metrics.centrality_z_scores:
            # Save detailed z-scores
            with open(self.output_dir / f"{prefix}_centrality_significance.json", 'w') as f:
                json.dump(metrics.centrality_z_scores, f, indent=2, default=str)
            
            # Save significant nodes CSV
            sig_nodes = metrics.centrality_z_scores.get('significant_nodes', {})
            if any(sig_nodes.values()):
                with open(self.output_dir / f"{prefix}_significant_nodes.csv", 'w') as f:
                    f.write("node,significance_level,betweenness_z,closeness_z,degree_z,eigenvector_z\n")
                    
                    all_z_scores = {
                        'betweenness_z': metrics.centrality_z_scores.get('betweenness_z', {}),
                        'closeness_z': metrics.centrality_z_scores.get('closeness_z', {}),
                        'degree_z': metrics.centrality_z_scores.get('degree_z', {}),
                        'eigenvector_z': metrics.centrality_z_scores.get('eigenvector_z', {})
                    }
                    
                    # Write highly significant nodes
                    for node in sig_nodes.get('highly_significant', []):
                        z_vals = [all_z_scores[measure].get(node, 0.0) for measure in 
                                ['betweenness_z', 'closeness_z', 'degree_z', 'eigenvector_z']]
                        f.write(f"{node},highly_significant,{','.join(f'{z:.3f}' for z in z_vals)}\n")
                    
                    # Write significant nodes
                    for node in sig_nodes.get('significant', []):
                        z_vals = [all_z_scores[measure].get(node, 0.0) for measure in 
                                ['betweenness_z', 'closeness_z', 'degree_z', 'eigenvector_z']]
                        f.write(f"{node},significant,{','.join(f'{z:.3f}' for z in z_vals)}\n")
            
            # Save significance summary
            with open(self.output_dir / f"{prefix}_significance_summary.txt", 'w') as f:
                f.write("Centrality Significance Analysis Summary\n")
                f.write("=" * 45 + "\n")
                
                summary = metrics.centrality_z_scores.get('summary', {})
                f.write(f"Highly significant nodes (|z| > 2.5): {summary.get('n_highly_significant', 0)}\n")
                f.write(f"Significant nodes (|z| > 1.96): {summary.get('n_significant', 0)}\n")
                f.write(f"Notable nodes (|z| > 1.0): {summary.get('n_notable', 0)}\n")
                f.write(f"Total statistically significant: {summary.get('total_significant', 0)}\n")
                
                # List top significant nodes
                highly_sig = sig_nodes.get('highly_significant', [])
                if highly_sig:
                    f.write(f"\nHighly Significant Nodes:\n")
                    for node in highly_sig[:10]:
                        f.write(f"  {node}\n")
            
            # Create advanced network visualization
            self._create_advanced_network_visualization(simulation, prefix)
    
    def _create_dccm_visualization(self, dccm_matrix: np.ndarray, prefix: str):
        """Create DCCM heatmap visualization"""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            plt.figure(figsize=(10, 8))
            
            # Create heatmap with diverging colormap
            sns.heatmap(dccm_matrix, 
                       cmap='RdBu_r', 
                       center=0,
                       vmin=-1, vmax=1,
                       square=True,
                       cbar_kws={'label': 'Cross-Correlation'})
            
            plt.title('Dynamic Cross-Correlation Matrix (DCCM)')
            plt.xlabel('Residue Index')
            plt.ylabel('Residue Index')
            plt.tight_layout()
            
            # Save plot
            plt.savefig(self.output_dir / f"{prefix}_dccm_heatmap.png", dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not create DCCM visualization: {e}")
    
    def _create_pca_visualization(self, pca_data: Dict[str, Any], prefix: str):
        """Create PCA analysis visualization"""
        try:
            import matplotlib.pyplot as plt
            
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
            
            # Plot 1: Eigenvalue spectrum (scree plot)
            eigenvals = pca_data['pca_eigenvalues']
            ax1.plot(range(1, len(eigenvals) + 1), eigenvals, 'bo-')
            ax1.set_xlabel('Principal Component')
            ax1.set_ylabel('Eigenvalue')
            ax1.set_title('PCA Eigenvalue Spectrum')
            ax1.grid(True, alpha=0.3)
            
            # Plot 2: Variance explained
            var_explained = pca_data['pca_variance_explained'] * 100
            ax2.bar(range(1, len(var_explained) + 1), var_explained)
            ax2.set_xlabel('Principal Component')
            ax2.set_ylabel('Variance Explained (%)')
            ax2.set_title('Variance Explained by Each PC')
            ax2.grid(True, alpha=0.3)
            
            # Plot 3: Cumulative variance explained
            cum_var = pca_data['pca_cumulative_variance'] * 100
            ax3.plot(range(1, len(cum_var) + 1), cum_var, 'ro-')
            ax3.axhline(y=80, color='k', linestyle='--', alpha=0.5, label='80%')
            ax3.set_xlabel('Principal Component')
            ax3.set_ylabel('Cumulative Variance Explained (%)')
            ax3.set_title('Cumulative Variance Explained')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
            
            # Plot 4: PC1 vs PC2 trajectory projection
            if 'principal_components' in pca_data and pca_data['principal_components'] is not None:
                pc_proj = pca_data['principal_components']
                if pc_proj.shape[1] >= 2:
                    scatter = ax4.scatter(pc_proj[:, 0], pc_proj[:, 1], 
                                        c=range(len(pc_proj[:, 0])), cmap='viridis', alpha=0.6)
                    ax4.set_xlabel(f'PC1 ({var_explained[0]:.1f}%)')
                    ax4.set_ylabel(f'PC2 ({var_explained[1]:.1f}%)')
                    ax4.set_title('Trajectory Projection on PC1-PC2')
                    plt.colorbar(scatter, ax=ax4, label='Frame')
                else:
                    ax4.text(0.5, 0.5, 'Insufficient PCs\nfor projection plot', 
                           ha='center', va='center', transform=ax4.transAxes)
            else:
                ax4.text(0.5, 0.5, 'Principal components\nnot available', 
                        ha='center', va='center', transform=ax4.transAxes)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / f"{prefix}_pca_analysis.png", dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not create PCA visualization: {e}")
    
    def _create_landscape_visualization(self, landscape_data: Dict[str, Any], prefix: str):
        """Create comprehensive energy landscape visualizations"""
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            from matplotlib.colors import LinearSegmentedColormap
            
            landscape = landscape_data['energy_landscape']
            pc1_bins = landscape_data['landscape_pc1_bins']
            pc2_bins = landscape_data['landscape_pc2_bins']
            gradient_x = landscape_data['landscape_gradient_x']
            gradient_y = landscape_data['landscape_gradient_y']
            laplacian = landscape_data['landscape_laplacian']
            minima = landscape_data.get('landscape_minima', [])
            
            # Create comprehensive landscape visualization
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
            
            # Plot 1: Energy landscape with contours
            pc1_centers = (pc1_bins[:-1] + pc1_bins[1:]) / 2
            pc2_centers = (pc2_bins[:-1] + pc2_bins[1:]) / 2
            PC1, PC2 = np.meshgrid(pc1_centers, pc2_centers)
            
            # Energy contour plot
            contour_levels = np.linspace(0, min(30, np.max(landscape)), 15)
            cs = ax1.contourf(PC1, PC2, landscape, levels=contour_levels, cmap='viridis')
            contour_lines = ax1.contour(PC1, PC2, landscape, levels=contour_levels, 
                                      colors='white', alpha=0.3, linewidths=0.5)
            ax1.clabel(contour_lines, inline=True, fontsize=8, fmt='%.1f')
            
            # Mark minima
            if minima:
                for i, (pc1_min, pc2_min, energy_min) in enumerate(minima[:5]):
                    ax1.plot(pc1_min, pc2_min, 'r*', markersize=15, markeredgewidth=1, markeredgecolor='white')
                    ax1.text(pc1_min, pc2_min + 0.1, f'M{i+1}', ha='center', va='bottom', 
                           color='white', fontweight='bold', fontsize=9)
            
            ax1.set_xlabel('PC1')
            ax1.set_ylabel('PC2')
            ax1.set_title('Free Energy Landscape (kJ/mol)')
            
            # Add colorbar
            cbar1 = plt.colorbar(cs, ax=ax1, shrink=0.8)
            cbar1.set_label('Energy (kJ/mol)')
            
            # Plot 2: 3D surface plot
            from mpl_toolkits.mplot3d import Axes3D
            ax2.remove()
            ax2 = fig.add_subplot(2, 2, 2, projection='3d')
            
            # Downsample for 3D plot if too dense
            step = max(1, landscape.shape[0] // 25)
            X_sub = PC1[::step, ::step]
            Y_sub = PC2[::step, ::step]
            Z_sub = landscape[::step, ::step]
            
            surf = ax2.plot_surface(X_sub, Y_sub, Z_sub, cmap='viridis', alpha=0.8, 
                                  linewidth=0, antialiased=True)
            ax2.set_xlabel('PC1')
            ax2.set_ylabel('PC2')
            ax2.set_zlabel('Energy (kJ/mol)')
            ax2.set_title('3D Energy Surface')
            ax2.view_init(elev=30, azim=45)
            
            # Plot 3: Gradient magnitude (force field)
            gradient_mag = np.sqrt(gradient_x**2 + gradient_y**2)
            im3 = ax3.imshow(gradient_mag, extent=[pc1_bins[0], pc1_bins[-1], pc2_bins[0], pc2_bins[-1]], 
                           origin='lower', cmap='plasma', aspect='auto')
            
            # Add gradient vectors (downsampled)
            step = max(1, landscape.shape[0] // 15)
            X_grad = PC1[::step, ::step]
            Y_grad = PC2[::step, ::step]
            U_grad = -gradient_x[::step, ::step]  # Negative for forces (point downhill)
            V_grad = -gradient_y[::step, ::step]
            
            ax3.quiver(X_grad, Y_grad, U_grad, V_grad, scale=None, alpha=0.7, color='white', width=0.002)
            ax3.set_xlabel('PC1')
            ax3.set_ylabel('PC2')
            ax3.set_title('Gradient Magnitude & Force Vectors')
            
            cbar3 = plt.colorbar(im3, ax=ax3, shrink=0.8)
            cbar3.set_label('|∇G| (kJ/mol/unit)')
            
            # Plot 4: Laplacian (curvature)
            im4 = ax4.imshow(laplacian, extent=[pc1_bins[0], pc1_bins[-1], pc2_bins[0], pc2_bins[-1]], 
                           origin='lower', cmap='RdBu_r', aspect='auto')
            
            # Mark minima on Laplacian plot
            if minima:
                for i, (pc1_min, pc2_min, _) in enumerate(minima[:5]):
                    ax4.plot(pc1_min, pc2_min, 'ko', markersize=8, markerfacecolor='yellow', 
                           markeredgewidth=2, markeredgecolor='black')
            
            ax4.set_xlabel('PC1')
            ax4.set_ylabel('PC2')
            ax4.set_title('Laplacian (Curvature)')
            
            cbar4 = plt.colorbar(im4, ax=ax4, shrink=0.8)
            cbar4.set_label('∇²G (curvature)')
            
            plt.tight_layout()
            plt.savefig(self.output_dir / f"{prefix}_energy_landscape.png", dpi=300, bbox_inches='tight')
            plt.close()
            
            # Create separate detailed contour plot
            fig2, ax = plt.subplots(1, 1, figsize=(10, 8))
            
            # High-resolution contour plot
            contour_levels_detailed = np.linspace(0, min(25, np.max(landscape)), 25)
            cs_detailed = ax.contourf(PC1, PC2, landscape, levels=contour_levels_detailed, cmap='viridis')
            contour_lines_detailed = ax.contour(PC1, PC2, landscape, levels=contour_levels_detailed, 
                                             colors='white', alpha=0.4, linewidths=0.8)
            
            # Add energy minima with labels
            if minima:
                for i, (pc1_min, pc2_min, energy_min) in enumerate(minima):
                    ax.plot(pc1_min, pc2_min, 'r*', markersize=20, markeredgewidth=2, markeredgecolor='white')
                    ax.text(pc1_min, pc2_min + 0.15, f'Min {i+1}\n{energy_min:.1f} kJ/mol', 
                           ha='center', va='bottom', color='white', fontweight='bold', 
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.7))
            
            ax.set_xlabel('PC1 Projection', fontsize=12)
            ax.set_ylabel('PC2 Projection', fontsize=12)
            ax.set_title(f'Free Energy Landscape - Temperature: {landscape_data.get("landscape_temperature", 310):.0f} K', 
                        fontsize=14)
            
            # Enhanced colorbar
            cbar = plt.colorbar(cs_detailed, ax=ax, shrink=0.8, pad=0.02)
            cbar.set_label('Free Energy (kJ/mol)', fontsize=12)
            cbar.ax.tick_params(labelsize=10)
            
            # Add some statistics as text
            stats_text = f"Energy Range: {np.max(landscape) - np.min(landscape):.1f} kJ/mol\n"
            stats_text += f"Minima Found: {len(minima)}\n"
            stats_text += f"Frames: {landscape_data.get('landscape_total_frames', 'N/A')}"
            
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8),
                   fontsize=10)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / f"{prefix}_energy_landscape_detailed.png", dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not create energy landscape visualization: {e}")
            import traceback
            print(traceback.format_exc())
    
    def _create_advanced_network_visualization(self, simulation: MDSimulation, prefix: str):
        """Create advanced network analysis visualizations"""
        if not simulation.network_metrics:
            return
            
        metrics = simulation.network_metrics
        
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as patches
            
            # Create comprehensive network analysis figure
            fig = plt.figure(figsize=(20, 15))
            
            # 1. Community structure visualization
            ax1 = plt.subplot(3, 3, 1)
            if hasattr(metrics, 'communities_detailed') and metrics.communities_detailed:
                community_sizes = metrics.communities_detailed.get('community_sizes', [])
                if community_sizes:
                    ax1.bar(range(len(community_sizes)), community_sizes)
                    ax1.set_xlabel('Community ID')
                    ax1.set_ylabel('Community Size')
                    ax1.set_title('Community Size Distribution')
                    ax1.grid(True, alpha=0.3)
            
            # 2. Centrality significance
            ax2 = plt.subplot(3, 3, 2) 
            if hasattr(metrics, 'centrality_z_scores') and metrics.centrality_z_scores:
                summary = metrics.centrality_z_scores.get('summary', {})
                categories = ['Highly Sig.', 'Significant', 'Notable']
                counts = [summary.get('n_highly_significant', 0), 
                         summary.get('n_significant', 0),
                         summary.get('n_notable', 0)]
                
                colors = ['red', 'orange', 'yellow']
                ax2.bar(categories, counts, color=colors, alpha=0.7)
                ax2.set_ylabel('Number of Nodes')
                ax2.set_title('Centrality Significance')
                ax2.grid(True, alpha=0.3)
            
            # 3. Path length distribution
            ax3 = plt.subplot(3, 3, 3)
            if hasattr(metrics, 'path_metrics') and metrics.path_metrics:
                path_dist = metrics.path_metrics.get('path_length_distribution', {})
                if path_dist and 'lengths' in path_dist:
                    lengths = path_dist['lengths']
                    counts = path_dist['counts']
                    ax3.bar(lengths, counts, alpha=0.7)
                    ax3.set_xlabel('Path Length')
                    ax3.set_ylabel('Frequency')
                    ax3.set_title('Path Length Distribution')
                    ax3.grid(True, alpha=0.3)
            
            # 4. Network robustness
            ax4 = plt.subplot(3, 3, 4)
            if hasattr(metrics, 'network_robustness') and metrics.network_robustness:
                nr = metrics.network_robustness
                robustness_measures = ['Node Attack', 'Random Failure', 'Edge Attack']
                robustness_values = [
                    nr.get('node_attack_robustness', 0),
                    nr.get('random_failure_robustness', 0),
                    nr.get('edge_attack_robustness', 0)
                ]
                
                bars = ax4.bar(robustness_measures, robustness_values, 
                              color=['red', 'blue', 'green'], alpha=0.7)
                ax4.set_ylabel('Robustness Score')
                ax4.set_title('Network Robustness')
                ax4.set_ylim(0, 1)
                ax4.grid(True, alpha=0.3)
                
                # Add value labels on bars
                for bar, value in zip(bars, robustness_values):
                    height = bar.get_height()
                    ax4.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                            f'{value:.3f}', ha='center', va='bottom')
            
            # 5. Allosteric hotspots
            ax5 = plt.subplot(3, 3, 5)
            if hasattr(metrics, 'allosteric_pathways') and metrics.allosteric_pathways:
                hotspots = metrics.allosteric_pathways.get('allosteric_hotspots', [])
                if hotspots:
                    top_hotspots = hotspots[:10]  # Top 10
                    nodes = [h['node'] for h in top_hotspots]
                    scores = [h['hotspot_score'] for h in top_hotspots]
                    
                    ax5.barh(range(len(nodes)), scores)
                    ax5.set_yticks(range(len(nodes)))
                    ax5.set_yticklabels(nodes, fontsize=8)
                    ax5.set_xlabel('Hotspot Score')
                    ax5.set_title('Top Allosteric Hotspots')
                    ax5.grid(True, alpha=0.3)
            
            # 6. Communication efficiency heatmap
            ax6 = plt.subplot(3, 3, 6)
            if hasattr(metrics, 'allosteric_pathways') and metrics.allosteric_pathways:
                comm_eff = metrics.allosteric_pathways.get('communication_efficiency', {})
                if comm_eff and len(comm_eff) > 1:
                    # Create a small efficiency matrix for top pairs
                    top_pairs = sorted(comm_eff.items(), key=lambda x: x[1], reverse=True)[:25]
                    
                    # Extract unique nodes
                    all_nodes = set()
                    for pair_key, _ in top_pairs:
                        source, target = pair_key.split('_to_')
                        all_nodes.add(source)
                        all_nodes.add(target)
                    
                    nodes_list = sorted(list(all_nodes))[:10]  # Limit to 10 nodes
                    
                    if len(nodes_list) >= 2:
                        # Create efficiency matrix
                        eff_matrix = np.zeros((len(nodes_list), len(nodes_list)))
                        
                        for i, source in enumerate(nodes_list):
                            for j, target in enumerate(nodes_list):
                                key1 = f"{source}_to_{target}"
                                key2 = f"{target}_to_{source}"
                                eff = comm_eff.get(key1, comm_eff.get(key2, 0))
                                eff_matrix[i, j] = eff
                        
                        im = ax6.imshow(eff_matrix, cmap='viridis', aspect='auto')
                        ax6.set_xticks(range(len(nodes_list)))
                        ax6.set_yticks(range(len(nodes_list)))
                        ax6.set_xticklabels(nodes_list, rotation=45, fontsize=8)
                        ax6.set_yticklabels(nodes_list, fontsize=8)
                        ax6.set_title('Communication Efficiency')
                        plt.colorbar(im, ax=ax6, shrink=0.8)
            
            # 7. Centrality correlation plot
            ax7 = plt.subplot(3, 3, 7)
            if metrics.betweenness_centrality and metrics.degree_centrality:
                nodes = list(set(metrics.betweenness_centrality.keys()) & 
                           set(metrics.degree_centrality.keys()))
                if nodes:
                    bet_vals = [metrics.betweenness_centrality[node] for node in nodes]
                    deg_vals = [metrics.degree_centrality[node] for node in nodes]
                    
                    ax7.scatter(deg_vals, bet_vals, alpha=0.6)
                    ax7.set_xlabel('Degree Centrality')
                    ax7.set_ylabel('Betweenness Centrality')
                    ax7.set_title('Centrality Correlation')
                    ax7.grid(True, alpha=0.3)
                    
                    # Add correlation coefficient
                    if len(bet_vals) > 1:
                        corr_coef = np.corrcoef(deg_vals, bet_vals)[0, 1]
                        ax7.text(0.05, 0.95, f'r = {corr_coef:.3f}', 
                                transform=ax7.transAxes, fontsize=10,
                                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # 8. Network efficiency metrics
            ax8 = plt.subplot(3, 3, 8)
            if hasattr(metrics, 'path_metrics') and metrics.path_metrics:
                pm = metrics.path_metrics
                efficiency_types = ['Global', 'Local', 'Network']
                efficiency_values = [
                    pm.get('global_efficiency', 0),
                    pm.get('local_efficiency', 0), 
                    pm.get('network_efficiency', 0)
                ]
                
                bars = ax8.bar(efficiency_types, efficiency_values, 
                              color=['purple', 'cyan', 'orange'], alpha=0.7)
                ax8.set_ylabel('Efficiency')
                ax8.set_title('Network Efficiency Measures')
                ax8.grid(True, alpha=0.3)
                
                # Add value labels
                for bar, value in zip(bars, efficiency_values):
                    height = bar.get_height()
                    ax8.text(bar.get_x() + bar.get_width()/2., height + 0.001,
                            f'{value:.4f}', ha='center', va='bottom', fontsize=9)
            
            # 9. Summary statistics text
            ax9 = plt.subplot(3, 3, 9)
            ax9.axis('off')
            
            summary_text = f"Advanced Network Analysis Summary\n"
            summary_text += f"{'='*40}\n\n"
            
            if hasattr(metrics, 'communities_detailed') and metrics.communities_detailed:
                cd = metrics.communities_detailed
                summary_text += f"Communities: {cd.get('n_communities', 0)}\n"
                summary_text += f"Modularity: {cd.get('modularity', 0):.4f}\n\n"
            
            if hasattr(metrics, 'path_metrics') and metrics.path_metrics:
                pm = metrics.path_metrics
                summary_text += f"Path Length: {pm.get('characteristic_path_length', 0):.3f}\n"
                summary_text += f"Diameter: {pm.get('network_diameter', 0):.3f}\n\n"
            
            if hasattr(metrics, 'allosteric_pathways') and metrics.allosteric_pathways:
                ap = metrics.allosteric_pathways
                pathways = ap.get('pathways', [])
                hotspots = ap.get('allosteric_hotspots', [])
                summary_text += f"Allosteric Pathways: {len(pathways)}\n"
                summary_text += f"Hotspot Residues: {len(hotspots)}\n\n"
            
            if hasattr(metrics, 'network_robustness') and metrics.network_robustness:
                nr = metrics.network_robustness
                summary_text += f"Network Resilience: {nr.get('network_resilience', 0):.3f}\n"
                
                frag_point = nr.get('fragmentation_point')
                if frag_point:
                    summary_text += f"Fragmentation Point: {frag_point} nodes\n"
            
            ax9.text(0.05, 0.95, summary_text, transform=ax9.transAxes, 
                    fontsize=10, verticalalignment='top', fontfamily='monospace',
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig(self.output_dir / f"{prefix}_advanced_network_analysis.png", 
                       dpi=300, bbox_inches='tight')
            plt.close()
            
        except Exception as e:
            print(f"Warning: Could not create advanced network visualization: {e}")
            import traceback
            print(traceback.format_exc())


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
            
            # Compute dynamic analysis (DCCM and PCA)
            dynamic_results = self.analyzer.compute_dynamic_analysis(simulation)
            
            # Store dynamic results in simulation
            if dynamic_results:
                simulation.dynamic_analysis = dynamic_results
            
            # Trajectory segmentation analysis
            if self.analysis_config.segments > 1:
                segment_results = self.analyzer.analyze_trajectory_segments(simulation)
            
            # Save individual results
            self.output_manager.save_simulation_results(simulation)
            
            results_summary[name] = {
                'contact_maps_computed': list(contact_maps.keys()),
                'network_nodes': network_metrics.n_nodes,
                'network_edges': network_metrics.n_edges,
                'network_density': network_metrics.density,
                'dynamic_analysis_computed': bool(dynamic_results),
                'dccm_computed': 'dccm_matrix' in dynamic_results if dynamic_results else False,
                'pca_computed': 'pca_eigenvalues' in dynamic_results if dynamic_results else False
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
