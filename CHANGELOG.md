# Changelog

All notable changes to MD-Compare will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Integration with protein structure databases (PDB, AlphaFold)
- Machine learning models 
- Support for additional interaction types (π-π stacking, cation-π)
- Enhanced visualization 

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

## [0.9.0] - 2024-12-XX (Development/Pre-release)

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

## Migration Guide

### From Development to v1.0.0
- Adopt new class-based architecture (MDSimulation, NetworkAnalyzer, MDComparator)
- Update to new CLI command structure (`md-compare single/compare/diff`)
- Migrate to JSON-based configuration files for multi-simulation analyses
- Update import statements to use new module structure

### Configuration File Format Changes
- New JSON schema with `simulations` and `analysis` sections
- Enhanced metadata support for tracking simulation conditions
- Flexible parameter specification with reasonable defaults

---

**Note**: This represents the initial major release establishing MD-Compare as a 
comprehensive toolkit for molecular dynamics simulation comparison. Future 
versions will build upon this foundation with additional analysis methods, 
enhanced visualization capabilities, and integration with other computational 
biology tools.

