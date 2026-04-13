# Changelog

All notable changes to MD-Compare will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- pyEMMA integration 
- Pocket detection and characterization across simulations

## [1.3.0] - 2026-04-13

### Added
- **Advanced Community Detection**:
  - Multi-algorithm support: Leiden, Louvain, Spectral, and Hierarchical clustering
  - Modularity optimization with quality assessment
  - Intra/inter-community density calculations
  - Community size distributions and node assignments

- **Comprehensive Path Metrics Analysis**:
  - Global and local efficiency measurements
  - Characteristic path length and network diameter calculation
  - Path length distributions with statistical analysis
  - Node eccentricity computation (remoteness measures)
  - Critical pathway identification through network hubs

- **Sophisticated Allosteric Pathway Analysis**:
  - Automated functional region detection (HIV protease-specific)
  - Communication efficiency quantification between functional sites
  - Allosteric hotspot discovery with frequency-based scoring
  - Pathway bottleneck identification and redundancy assessment
  - Inter-community bridge analysis for functional regions

- **Network Robustness Analysis**:
  - Attack tolerance testing (targeted, random, edge-based)
  - Giant component tracking and fragmentation point identification
  - Critical node and edge identification for connectivity
  - Network resilience scoring and vulnerability assessment

- **Statistical Significance Analysis**:
  - Centrality z-score computation for all measures
  - Statistical outlier identification (highly significant, significant, notable)
  - Multi-measure analysis across betweenness, closeness, degree, eigenvector
  - Node classification by statistical significance levels

- **Advanced Visualization Dashboard**:
  - 9-panel comprehensive network analysis overview
  - Community structure, significance, and robustness visualizations
  - Allosteric hotspot rankings and communication efficiency heatmaps
  - Path analysis and centrality correlation plots
  - Summary statistics panel with key metrics

### New CLI Options
- `--community-method leiden|louvain|spectral|hierarchical` - Choose community detection algorithm
- `--no-communities` - Skip community detection
- `--no-paths` - Skip detailed path analysis
- `--no-allosteric` - Skip allosteric pathway analysis
- `--allosteric-sources A_50 B_50` - Specify allosteric source nodes
- `--allosteric-targets A_25 B_25` - Specify allosteric target nodes

### Enhanced Output Files
- `*_community_analysis.json` - Detailed community detection results
- `*_community_assignments.csv` - Node-to-community mapping
- `*_path_metrics.json` - Comprehensive path analysis data
- `*_allosteric_analysis.json` - Complete allosteric pathway analysis
- `*_allosteric_hotspots.csv` - Ranked hotspot residues with scores
- `*_robustness_analysis.json` - Network resilience data
- `*_centrality_significance.json` - Statistical significance analysis
- `*_significant_nodes.csv` - Statistically significant nodes
- `*_advanced_network_analysis.png` - Comprehensive visualization dashboard

### Scientific Enhancements
- **HIV Protease Specialization**: Automated detection of flap regions, active site, and interface
- **Allosteric Mechanism Mapping**: Direct identification of communication pathways
- **Drug Target Discovery**: Network-based identification of critical therapeutic targets
- **Mutation Impact Prediction**: Robustness analysis for resistance mechanism studies

### Performance Improvements
- **Efficient Algorithms**: Optimized community detection with multiple fallback methods
- **Scalable Analysis**: Smart sampling for large networks (>500 nodes)
- **Memory Optimization**: Efficient sparse matrix operations for path analysis
- **Progress Monitoring**: Detailed progress reporting for long computations

### Technical Architecture
- **Modular Design**: Separate methods for each analysis component
- **Error Handling**: Graceful fallback when advanced algorithms fail
- **Statistical Rigor**: Proper significance testing with z-score analysis
- **Cross-Platform**: Compatible with Windows, macOS, and Linux

### Dependencies
- **Optional Enhancement**: `igraph` + `leidenalg` for Leiden community detection
- **Scikit-learn**: For spectral clustering and advanced network analysis
- **Scipy Enhanced**: Additional clustering and statistical functions

### Usage Examples
```bash
# Full advanced network analysis
md-compare single -t system.pdb -x traj.xtc -n analysis \
  --community-method leiden \
  --allosteric-sources A_50 B_50 \
  --allosteric-targets A_25 B_25

# Fast analysis (skip advanced features)
md-compare single -t system.pdb -x traj.xtc -n analysis \
  --no-allosteric --no-paths --community-method louvain
```



## [1.2.0] - 2026-04-10

### Added
- **Dynamic Analysis Capabilities**: 
  - Dynamic Cross-Correlation Matrix (DCCM) analysis for residue motion correlations
  - Principal Component Analysis (PCA) for identifying dominant motion modes
  - Automatic integration with existing network analysis workflow

- **Energy Landscape Analysis**:
  - Free energy surface calculation from PC1/PC2 projections using Boltzmann relation
  - Gaussian smoothing and thermodynamic analysis with temperature control
  - Energy minima detection with automated separation filtering
  - Energy barrier analysis between conformational states
  - Gradient and Laplacian computation for force field and curvature analysis
  - 2D histogram-based probability distributions with configurable binning

- **New CLI Options**:
  - `--compute-dccm` / `--no-dccm` flags for DCCM computation control (default: enabled)
  - `--compute-pca` / `--no-pca` flags for PCA analysis control (default: enabled)
  - `--pca-components N` to specify number of principal components (default: 10)
  - `--dccm-selection "selection"` to specify atoms for dynamic analysis (default: "name CA")
  - `--compute-landscape` / `--no-landscape` flags for energy landscape control (default: enabled)
  - `--landscape-temp 310` to set temperature for energy calculations (default: 310K)
  - `--landscape-bins 50` to specify histogram resolution (default: 50)
  - `--landscape-sigma 1.0` to control Gaussian smoothing (default: 1.0)

- **Enhanced Output Files**:
  - `*_dccm_matrix.npy/.csv` - Dynamic cross-correlation matrices
  - `*_dccm_heatmap.png` - DCCM visualization heatmaps
  - `*_pca_eigenvalues.npy` - PCA eigenvalues and variance explained
  - `*_pca_eigenvectors.npy` - PCA eigenvectors (motion modes)
  - `*_pca_variance_explained.npy` - Variance explained by each component
  - `*_principal_components.npy` - Trajectory projections onto PCs
  - `*_pca_analysis.png` - Comprehensive PCA plots (scree, variance, projections)
  - `*_energy_landscape.npy` - Free energy surface matrix
  - `*_landscape_pc1_bins.npy/.._pc2_bins.npy` - Histogram bin edges
  - `*_landscape_gradient_x/y.npy` - Energy gradients (force field components)
  - `*_landscape_laplacian.npy` - Curvature analysis matrix
  - `*_energy_landscape.png` - 4-panel comprehensive landscape analysis
  - `*_energy_landscape_detailed.png` - High-resolution contour plot with statistics
  - `*_landscape_summary.txt` - Energy landscape analysis summary with minima/barriers
  - `*_pca_summary.txt` - PCA analysis summary with variance breakdown
  - `*_dynamic_summary.txt` - Overall dynamic analysis summary

- **Advanced Visualization**:
  - Automatic DCCM heatmap generation with diverging colormap
  - Multi-panel PCA analysis plots (eigenvalue spectrum, variance explained, trajectory projections)
  - Color-coded correlation matrices showing positive/negative correlations
  - 4-panel energy landscape visualization: contour plot, 3D surface, gradient field, curvature
  - High-resolution detailed contour plots with energy minima marked
  - PC1-PC2 trajectory evolution plots with temporal coloring

- **Scientific Energy Landscape Features**:
  - Boltzmann relation implementation: G = -kB*T*ln(P) for thermodynamically correct energies
  - Automated energy minima detection with configurable separation filtering
  - Energy barrier height analysis between conformational states
  - Force field computation from energy gradients for pathway analysis
  - Laplacian curvature analysis for stability/instability region identification
  - Temperature-dependent conformational accessibility assessment
  
### Changed
- **AnalysisConfig**: Extended with dynamic analysis configuration options
  - Added `compute_dccm`, `compute_pca`, `pca_components`, `dccm_selection` parameters
  - Added `landscape_temperature`, `landscape_bins`, `landscape_sigma` controls
- **NetworkMetrics**: Extended dataclass to include DCCM, PCA, and energy landscape results
- **CLI Interface**: All analysis modes (single, compare, diff) now support dynamic analysis flags
- **Output Management**: Automatic saving and visualization of dynamic analysis results

### Enhanced
- **Example 1**: Updated single simulation example with DCCM and PCA usage demonstrations
- **Documentation**: Added dynamic analysis interpretation guidelines and advanced analysis scripts
- **Error Handling**: Graceful fallback when dynamic analysis fails due to insufficient data
- **Memory Management**: Efficient coordinate collection and matrix operations for large systems

### Performance
- **Coordinate Collection**: Optimized trajectory iteration for dynamic analysis
- **Matrix Operations**: Efficient DCCM computation with vectorized operations
- **Progress Monitoring**: Detailed progress reporting for long DCCM computations
- **Memory Optimization**: Smart handling of large coordinate arrays

### Technical Improvements
- **PCA Implementation**: Robust eigenvalue decomposition with  DCCM, PCA, and energy landscape demonstrations
- **DCCM Algorithm**: Mathematically correct cross-correlation computation
- **Preprocessing Integration**: Proper alignment and centering for dynamic analysis
- **Energy Landscape**: Scientifically validated implementation following established methodologies
- **Visualization Pipeline**: Automated generation of publication-quality plots

### Dependencies
- **SciPy**: Enhanced requirement for robust linear algebra operations (`scipy.linalg.eigh`)
- **Matplotlib/Seaborn**: Leveraged for advanced visualization capabilities
- **NumPy**: Optimized matrix operations for efficient DCCM/PCA computation
- **Gaussian Filtering**: scipy.ndimage for energy landscape smoothing

### Usage Examples
```bash
# Enable all dynamic analysis (default behavior)
md-compare single -t system.pdb -x traj.xtc -n analysis

# Customize energy landscape parameters
md-compare single -t system.pdb -x traj.xtc -n analysis \
  --landscape-temp 300 --landscape-bins 60 --landscape-sigma 1.2

# High-resolution dynamics analysis
md-compare single -t system.pdb -x traj.xtc -n analysis \
  --pca-components 20 --dccm-selection "backbone" --landscape-bins 80

# Speed-optimized analysis (disable dynamics)
md-compare single -t system.pdb -x traj.xtc -n analysis \
  --no-dccm --no-pca --no-landscape
```

## [1.0.0] - 2026-04-08

### Added
- Initial release of MD-Compare: Molecular Dynamics Simulation Comparison Toolkit
- **Core Analysis Framework**:
  - MDSimulation class for individual simulation handling
  - NetworkAnalyzer class for contact maps and network metrics
  - MDComparator class for cross-simulation comparison
  - OutputManager class for comprehensive result generation
  - MDCompare class for workflow orchestration

- **Multi-Chain Protein Support**:
  - Automatic chain detection and proper residue mapping
  - Multi-chain contact map generation and visualization
  - Inter-chain vs intra-chain contact analysis
  - Chain boundary visualization in contact maps

- **Robust Network Analysis**:
  - Timeout protection for computationally expensive operations
  - Automatic algorithm switching for large networks (>500 nodes)
  - Multiple centrality measures: betweenness, closeness, eigenvector, degree
  - Community detection with modularity optimization
  - Path analysis with disconnected network handling
  - Network robustness and assortativity analysis

- **MD Preprocessing Capabilities**:
  - Trajectory alignment using user-defined selections
  - System centering to remove translational motion
  - Rigid body motion removal for network analysis
  - Support for custom alignment and centering selections

- **Statistical Robustness**:
  - Trajectory segmentation for confidence estimation
  - Contact persistence analysis and statistics
  - Error estimation through segment comparison
  - Performance monitoring and timing analysis

- **Comprehensive Comparison Features**:
  - Multi-simulation network property comparison
  - Residue-specific centrality comparisons
  - Contact pattern differential analysis
  - Statistical significance testing for network differences
  - Highly persistent contact identification

- **Multiple MD Format Support**:
  - GROMACS (.gro/.xtc/.trr) 
  - AMBER (.prmtop/.nc/.dcd)
  - CHARMM (.psf/.dcd)
  - NAMD (.psf/.dcd/.xsc)
  - Desmond (.cms/.dms/.dtr)
  - Generic PDB + trajectory combinations

- **Command-Line Interface**:
  - `md-compare single` - Single simulation analysis
  - `md-compare compare` - Multi-simulation comparison from JSON config
  - `md-compare diff` - Differential analysis between two simulations
  - `md-compare example-config` - Generate template configuration files

- **Advanced Output Generation**:
  - Contact frequency matrices (NumPy binary and CSV formats)
  - Network centrality measures with residue mapping
  - NetworkX-compatible graph files for external visualization
  - Comprehensive analysis reports in JSON format
  - Publication-ready contact map visualizations
  - Performance timing and memory usage reports

- **Interaction Type Support**:
  - Distance-based contacts with customizable cutoffs
  - Hydrogen bond detection and analysis
  - Salt bridge identification and quantification
  - Framework for additional interaction types

### Core Scripts and Modules
- `md_compare_core.py` - Main analysis framework with all core classes
- `md_compare_cli.py` - Command-line interface and workflow management  
- `utils.py` - Utility functions, preprocessing, and helper classes

### Configuration System
- JSON-based configuration for multi-simulation workflows
- Flexible analysis parameter specification
- Simulation metadata and description support
- Example configuration generators

### Documentation and Development
- `README.md` - Comprehensive usage guide and scientific applications
- `CONTRIBUTING.md` - Development guidelines and contribution processes
- `LICENSE` - MIT license with proper attribution requirements
- GitHub Actions CI/CD pipeline for automated testing
- Code quality enforcement (Black, flake8)
- Unit test framework setup

### Dependencies
- MDAnalysis >= 2.4.0 (molecular dynamics analysis)
- NetworkX >= 2.8 (network analysis algorithms)
- NumPy >= 1.21.0 (numerical computations)
- SciPy >= 1.7.0 (statistical analysis)
- Matplotlib >= 3.5.0 (visualization)
- Seaborn >= 0.11.0 (statistical plotting)
- Pandas >= 1.3.0 (data manipulation)

### Performance Characteristics
- Typical protein systems (100-500 residues): 1-15 minutes analysis time
- Large networks (>500 nodes): Automatic optimization with sampling
- Memory usage: <8GB for typical protein systems
- Support for 1000+ frame trajectories with segmentation

### Scientific Applications Validated
- Protein conformational change analysis
- Allosteric mechanism elucidation  
- Drug resistance mutation studies
- Multi-condition comparative studies
- Community structure and modular organization analysis

## [0.9.0] - 2026-04-06 (Development/Pre-release)

### Added
- Basic residue interaction network analysis framework
- Single-simulation contact map calculation
- NetworkX integration for basic centrality metrics
- Multi-chain protein system support foundation

### Issues Resolved in v1.0.0
- Comprehensive multi-simulation comparison capabilities
- Timeout protection for large network computations
- Statistical robustness through trajectory segmentation
- Performance optimization for large protein systems
- Advanced visualization and output generation
- Configuration-driven workflow management

---

## Version Numbering Scheme

- **Major versions (X.0.0)**: Breaking API changes, new core architecture
- **Minor versions (X.Y.0)**: New analysis methods, additional features
- **Patch versions (X.Y.Z)**: Bug fixes, performance improvements, documentation

## Support Policy

- **Current version**: Full support with active development
- **Previous major version**: Security updates and critical bug fixes  
- **Older versions**: Community support only

## Migration Notes

### From v1.0.0 to v1.1.0
- **Dynamic analysis enabled by default**: New installations will automatically compute DCCM and PCA
- **New output files**: Expect additional files in results directories (*_dccm_*, *_pca_*)
- **CLI compatibility**: All existing commands work unchanged; new flags are optional
- **Performance impact**: Analysis takes ~20-30% longer due to dynamic analysis
- **Memory requirements**: Increased memory usage for coordinate storage during DCCM/PCA

### For Existing Users:
- **Backward compatibility**: All v1.0.0 commands work identically in v1.1.0
- **Optional features**: Use `--no-dccm --no-pca` to replicate v1.0.0 behavior exactly
- **Output structure**: New files added to existing output directories without conflicts
- **Configuration files**: JSON configs automatically gain dynamic analysis options

### For Developers:
- **NetworkMetrics dataclass**: Extended with optional DCCM/PCA fields
- **MDSimulation objects**: New `dynamic_analysis` attribute available
- **Analysis workflow**: `compute_dynamic_analysis()` method available in NetworkAnalyzer
- **Output management**: Enhanced with dynamic analysis visualization methods

### Performance Optimization Tips:
```bash
# For large systems or long trajectories
md-compare single ... --dccm-selection "name CA" --pca-components 5

# For speed-critical analysis
md-compare single ... --no-dccm --no-pca

# For detailed dynamics analysis
md-compare single ... --pca-components 20 --dccm-selection "backbone"
```

---

**Note**: This represents the initial major release establishing MD-Compare as a 
comprehensive toolkit for molecular dynamics simulation comparison. Future 
versions will build upon this foundation with additional analysis methods, 
enhanced visualization capabilities, and integration with other computational 
biology tools.

