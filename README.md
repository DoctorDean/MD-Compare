# MD-Compare: Molecular Dynamics Simulation Comparison Toolkit

A comprehensive toolkit for comparing molecular dynamics simulations using network analysis, with a focus on understanding conformational changes and allosteric mechanisms in protein systems.

## 🧬 Overview

MD-Compare provides advanced computational tools for analyzing and comparing molecular dynamics simulations through residue interaction networks (RINs). The toolkit enables researchers to understand structural differences, identify key communication pathways, and quantify conformational changes between different simulation conditions.

### Key Capabilities

- **Multi-System Analysis**: Proper handling of multi-chain protein systems with automatic chain detection
- **Robust Network Metrics**: Centrality measures, community detection, and pathway analysis with timeout protection
- **MD Preprocessing**: Alignment and centering to remove rigid body motions
- **Comparative Analysis**: Direct comparison of multiple simulations with statistical significance testing
- **Statistical Robustness**: Trajectory segmentation for error estimation and confidence intervals
- **Multiple MD Formats**: Support for GROMACS, AMBER, CHARMM, NAMD, and Desmond trajectories

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/DoctorDean/md-compare.git
cd md-compare

# Install dependencies
pip install -r requirements.txt

# Install package
pip install -e .

# Verify installation
md-compare --help
```

### Basic Usage

```bash
# Single simulation analysis
md-compare single -t system.pdb -x trajectory.xtc -n my_simulation -o results/

# Compare multiple simulations
md-compare compare -c simulations_config.json -o comparison_results/

# Differential analysis between two simulations
md-compare diff -t1 wt.pdb -x1 wt.xtc -n1 wildtype \
                -t2 mut.pdb -x2 mut.xtc -n2 mutant -o diff_results/
```

## 📁 Repository Structure

```
├── md_compare_core.py              # Core analysis classes and workflow
├── md_compare_cli.py               # Command-line interface
├── utils.py                        # Utility functions and helpers
├── requirements.txt                # Python dependencies
├── setup.py                        # Package installation
├── LICENSE                         # MIT License
├── README.md                       # This file
├── examples/                       # Example configurations and scripts
├── docs/                           # Documentation
└── tests/                          # Unit tests
```

## 🔬 Scientific Applications

### Protein Conformational Analysis

```bash
# Analyze conformational changes in protein systems
md-compare single -t protein.pdb -x dynamics.xtc -n protein_analysis \
  --segments 5 --threshold 0.2
```

### Drug Resistance Mechanism Studies

```bash
# Compare wild-type vs mutant systems
md-compare diff -t1 wildtype.pdb -x1 wt_traj.xtc -n1 WT \
                -t2 mutant.pdb -x2 mut_traj.xtc -n2 MUT \
                --diff-threshold 0.1
```

### Allosteric Pathway Analysis

The toolkit automatically identifies:
- Communication pathways between distant sites
- Changes in network topology upon mutation
- Key residues for allosteric transmission
- Community structure and modular organization

### Multi-Condition Comparisons

```json
{
  "simulations": [
    {
      "name": "condition_A",
      "topology": "system_A.pdb",
      "trajectory": "traj_A.xtc",
      "description": "Control condition"
    },
    {
      "name": "condition_B", 
      "topology": "system_B.pdb",
      "trajectory": "traj_B.xtc",
      "description": "Treatment condition"
    }
  ]
}
```

## 📊 Output Files

### Network Analysis
- `*_centrality.csv` - Centrality measures for each residue
- `*_network.graphml` - Network files for visualization in Cytoscape/Gephi
- `*_system_info.txt` - System composition and analysis parameters

### Contact Analysis  
- `*_distance_contacts.npy` - Contact frequency matrices
- `*_contact_map.png` - Visual contact maps with chain boundaries
- `*_contacts.csv` - Contact data in tabular format

### Comparative Analysis
- `comparison_report.json` - Comprehensive comparison results
- `differential_*.json` - Detailed differential analysis
- `workflow_summary.json` - Complete analysis summary

## 🛠 Advanced Features

### Trajectory Segmentation

```bash
# Analyze trajectory in segments for statistical robustness
md-compare single -t system.pdb -x trajectory.xtc --segments 10 -o robust_analysis
```

### Custom Contact Definitions

```bash
# Multiple interaction types
md-compare single -t system.pdb -x trajectory.xtc \
  --interaction-types distance hbond salt_bridge \
  --cutoff 4.5 --threshold 0.15
```

### MD Preprocessing

```bash
# Disable preprocessing for pre-aligned trajectories
md-compare single -t system.pdb -x trajectory.xtc --no-preprocess

# Custom alignment selections
md-compare single -t system.pdb -x trajectory.xtc \
  --align-selection "name CA and resid 1-100" \
  --center-selection "protein"
```

### Performance Optimization

```bash
# Adjust timeout for large systems
md-compare single -t large_system.pdb -x trajectory.xtc \
  --timeout 600 --threshold 0.3
```

## 🧪 Configuration Files

Generate example configuration:
```bash
md-compare example-config -o my_simulations.json
```

Example configuration structure:
```json
{
  "simulations": [
    {
      "name": "simulation_1",
      "topology": "topology_1.pdb",
      "trajectory": "trajectory_1.xtc", 
      "selection": "protein and not name H*",
      "description": "Description of simulation 1",
      "metadata": {
        "temperature": 300,
        "force_field": "CHARMM36"
      }
    }
  ],
  "analysis": {
    "cutoffs": {"all_atom": 4.5, "ca_only": 8.0},
    "interaction_types": ["distance", "hbond"],
    "threshold": 0.2,
    "segments": 5,
    "preprocess": true
  }
}
```

## 📖 Documentation

### Key Parameters

- `--cutoff`: Distance cutoff for contacts in Å (default: 4.5)
- `--threshold`: Contact frequency cutoff for network edges (default: 0.2)  
- `--timeout`: Maximum time for expensive computations (default: 300s)
- `--segments`: Number of trajectory segments for statistics (default: 5)
- `--interaction-types`: Types of contacts to analyze (distance, hbond, salt_bridge)

### Troubleshooting

**Large network timeout issues:**
- Increase `--threshold` to reduce network density
- Decrease `--timeout` for faster approximate results
- Use distance-only analysis for initial exploration

**Multi-chain detection problems:**
- Check chain naming conventions in topology files
- Verify protein selection captures all desired chains
- Examine system info output for chain detection results

**Memory issues with large systems:**
- Use higher thresholds to create sparser networks
- Process shorter trajectory segments
- Consider distance-only analysis initially

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Development Setup

```bash
git clone https://github.com/DoctorDean/md-compare.git
cd md-compare
pip install -r requirements.txt
pip install -e .
python -m pytest tests/
```

## 📚 Citation

If you use MD-Compare in your research, please cite:

```bibtex
@software{md_compare_2026,
  title = {MD-Compare: Molecular Dynamics Simulation Comparison Toolkit},
  author = {Dean Sherry},
  year = {2026},
  url = {https://github.com/DoctorDean/md-compare}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **MDAnalysis** for trajectory analysis capabilities
- **NetworkX** for network analysis algorithms  
- **NumPy/SciPy** for numerical computations
- **Matplotlib** for visualization


**Disclaimer**: This software is for research purposes. Validate computational predictions through experimental methods when applicable.
