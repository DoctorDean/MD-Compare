# MD-Compare v1.4.0 Testing Guide

## **Post-Fix Validation Testing**

After applying the critical bug fixes for communication efficiency and organized output directories, follow this testing guide to ensure everything works correctly.

## ✅ **Quick Verification Test**

### **1. Basic Functionality Test**
```bash
# Test the fixed implementation
python md_compare_cli.py single \
  -t your_protein.pdb \
  -x your_trajectory.xtc \
  -n test_fixed_analysis \
  --compute-msm \
  --allosteric-sources A_50 B_50 \
  --allosteric-targets A_25 B_25
```

### **2. Expected Success Output**
```
Loading trajectory: your_trajectory.xtc
Topology: your_protein.pdb
Computing residue interaction network...
Network analysis complete: 198 nodes, 1010 edges
Centrality significance analysis complete: 32/198 statistically significant nodes
Computing PyEMMA Markov State Model analysis...
  MSM uses 96 of 100 total states
  Computing kinetic properties...
  ✓ Saved transition matrix as CSV: test_fixed_analysis_transition_matrix.csv
  ✓ Saved timescales as CSV: test_fixed_analysis_implied_timescales.csv
  ✓ Saved state populations as CSV: test_fixed_analysis_state_populations.csv
  Computing full communication efficiency matrix...
  ✓ Saved full communication efficiency matrix: test_fixed_analysis_communication_efficiency_matrix.csv
  ✓ Saved high-efficiency pairs: test_fixed_analysis_high_efficiency_pairs.csv
Output directory structure created at: test_fixed_analysis_results
✓ Analysis complete successfully
```

## **Directory Structure Validation**

After successful analysis, verify the organized output structure:

```bash
# Check directory structure (Windows)
dir test_fixed_analysis_results /s

# Check directory structure (Linux/macOS)
find test_fixed_analysis_results -type f | sort
```

### **Expected Structure**
```
test_fixed_analysis_results/
├── 01_residue_interaction_networks/
│   ├── test_fixed_analysis_advanced_network_analysis.png
│   ├── test_fixed_analysis_distance_contacts.csv
│   ├── test_fixed_analysis_network.graphml
│   ├── test_fixed_analysis_centrality.csv
│   ├── test_fixed_analysis_community_analysis.json
│   ├── test_fixed_analysis_community_assignments.csv
│   ├── test_fixed_analysis_robustness_analysis.json
│   └── test_fixed_analysis_significant_nodes.csv
├── 02_dynamic_cross_correlations/
│   ├── test_fixed_analysis_dccm_heatmap.png
│   ├── test_fixed_analysis_dccm_matrix.csv
│   └── test_fixed_analysis_dccm_matrix.npy
├── 03_principal_component_analysis/
│   ├── test_fixed_analysis_pca_analysis.png
│   ├── test_fixed_analysis_pca_eigenvalues.npy
│   ├── test_fixed_analysis_pca_eigenvectors.npy
│   └── test_fixed_analysis_pca_summary.txt
├── 04_energy_landscapes/
│   ├── test_fixed_analysis_energy_landscape.png       
│   ├── test_fixed_analysis_free_energy_surface.npy
│   └── test_fixed_analysis_landscape_summary.txt
├── 05_markov_state_models/
│   ├── test_fixed_analysis_msm_analysis.png
│   ├── test_fixed_analysis_transition_matrix.csv      # Excel-ready
│   ├── test_fixed_analysis_transition_matrix_summary.csv
│   ├── test_fixed_analysis_implied_timescales.csv
│   ├── test_fixed_analysis_state_populations.csv
│   ├── test_fixed_analysis_metastable_assignments.csv
│   ├── test_fixed_analysis_discrete_trajectory.npy
│   └── test_fixed_analysis_msm_summary.txt
├── 06_allosteric_pathways/                            
│   ├── test_fixed_analysis_communication_efficiency_matrix.csv  # ALL RESIDUE PAIRS
│   ├── test_fixed_analysis_high_efficiency_pairs.csv            # HIGH-EFFICIENCY SUMMARY
│   ├── test_fixed_analysis_allosteric_hotspots.csv
│   ├── test_fixed_analysis_allosteric_analysis.json
│   └── test_fixed_analysis_allosteric_summary.txt
└── 07_summary_reports/
    ├── test_fixed_analysis_complete_analysis_report.json
    ├── test_fixed_analysis_analysis_summary.txt
    ├── test_fixed_analysis_dynamic_summary.txt
    └── test_fixed_analysis_system_info.txt
```

## **Communication Efficiency Validation**

### **1. Check Full Matrix Size**
```python
import pandas as pd
import numpy as np

# Load the full communication efficiency matrix
matrix = pd.read_csv('test_fixed_analysis_results/06_allosteric_pathways/test_fixed_analysis_communication_efficiency_matrix.csv', index_col=0)

print(f"Communication efficiency matrix shape: {matrix.shape}")
print(f"Expected shape for 198 residues: (198, 198)")
print(f"Matrix covers all residue pairs: {matrix.shape[0] == 198 and matrix.shape[1] == 198}")

# Check for data completeness
non_zero_entries = (matrix > 0).sum().sum()
diagonal_entries = (np.diag(matrix) == 1.0).sum()
print(f"Non-zero entries: {non_zero_entries}")
print(f"Perfect self-communication (diagonal=1.0): {diagonal_entries}")
```

### **2. Validate High-Efficiency Pairs**
```python
# Load high-efficiency pairs
high_eff = pd.read_csv('test_fixed_analysis_results/06_allosteric_pathways/test_fixed_analysis_high_efficiency_pairs.csv')

print(f"Number of high-efficiency pairs (>0.2): {len(high_eff)}")
print("Top 5 most efficient communication pairs:")
print(high_eff.head())

# Check for cross-chain communication
cross_chain_pairs = high_eff[
    (high_eff['Source_Residue'].str.startswith('A_')) & 
    (high_eff['Target_Residue'].str.startswith('B_')) |
    (high_eff['Source_Residue'].str.startswith('B_')) & 
    (high_eff['Target_Residue'].str.startswith('A_'))
]
print(f"Cross-chain communication pairs: {len(cross_chain_pairs)}")
```

## **MSM Analysis Validation**

### **1. Transition Matrix Verification**
```python
# Load transition matrix
transition_matrix = pd.read_csv('test_fixed_analysis_results/05_markov_state_models/test_fixed_analysis_transition_matrix.csv', index_col=0)

print(f"Transition matrix shape: {transition_matrix.shape}")
print(f"Row sums (should be ~1.0): {transition_matrix.sum(axis=1).describe()}")
print(f"Matrix is stochastic: {np.allclose(transition_matrix.sum(axis=1), 1.0)}")

# Check connectivity
non_zero_transitions = (transition_matrix > 0).sum().sum()
total_possible = transition_matrix.shape[0] * transition_matrix.shape[1]
connectivity = non_zero_transitions / total_possible
print(f"Network connectivity: {connectivity:.3f}")
```

### **2. Kinetic Analysis Verification**
```python
# Load timescales
timescales = pd.read_csv('test_fixed_analysis_results/05_markov_state_models/test_fixed_analysis_implied_timescales.csv')

print("Implied timescales (frames):")
print(timescales.head(10))

# Load state populations  
populations = pd.read_csv('test_fixed_analysis_results/05_markov_state_models/test_fixed_analysis_state_populations.csv')

print(f"Most populated states:")
print(populations.head())
print(f"Population distribution sum: {populations['Population'].sum():.6f}")
```

## **Error Checking**

### **1. Check for Missing Files**
```bash
# Windows PowerShell
$expected_files = @(
    "06_allosteric_pathways/test_fixed_analysis_communication_efficiency_matrix.csv",
    "06_allosteric_pathways/test_fixed_analysis_high_efficiency_pairs.csv",
    "05_markov_state_models/test_fixed_analysis_transition_matrix.csv",
    "04_energy_landscapes/test_fixed_analysis_energy_landscape.png"
)

foreach ($file in $expected_files) {
    $path = "test_fixed_analysis_results/$file"
    if (Test-Path $path) {
        Write-Host "✓ Found: $file" -ForegroundColor Green
    } else {
        Write-Host "✗ Missing: $file" -ForegroundColor Red
    }
}
```

### **2. Validate File Sizes**
```python
import os

# Check critical files exist and have reasonable sizes
critical_files = {
    '06_allosteric_pathways/test_fixed_analysis_communication_efficiency_matrix.csv': (100000, 1000000),  # 100KB - 1MB
    '05_markov_state_models/test_fixed_analysis_transition_matrix.csv': (50000, 500000),   # 50KB - 500KB
    '06_allosteric_pathways/test_fixed_analysis_high_efficiency_pairs.csv': (1000, 50000), # 1KB - 50KB
}

for filepath, (min_size, max_size) in critical_files.items():
    full_path = f'test_fixed_analysis_results/{filepath}'
    if os.path.exists(full_path):
        size = os.path.getsize(full_path)
        if min_size <= size <= max_size:
            print(f"✓ {filepath}: {size:,} bytes (OK)")
        else:
            print(f"⚠ {filepath}: {size:,} bytes (unexpected size)")
    else:
        print(f"✗ {filepath}: File missing")
```

## 🧬 **Scientific Validation**

### **1. HIV Protease Communication Analysis**
```python
# Analyze cross-chain flap communication
high_eff = pd.read_csv('test_fixed_analysis_results/06_allosteric_pathways/test_fixed_analysis_high_efficiency_pairs.csv')

# HIV protease flap regions (typical residues)
flap_a = [f'A_{i}' for i in range(45, 55)]  # Chain A flap
flap_b = [f'B_{i}' for i in range(45, 55)]  # Chain B flap  
active_site = [f'A_{i}' for i in [25, 26, 27]] + [f'B_{i}' for i in [25, 26, 27]]

# Find flap-to-flap communication
flap_comm = high_eff[
    ((high_eff['Source_Residue'].isin(flap_a)) & (high_eff['Target_Residue'].isin(flap_b))) |
    ((high_eff['Source_Residue'].isin(flap_b)) & (high_eff['Target_Residue'].isin(flap_a)))
]

print(f"Flap-to-flap communication pathways: {len(flap_comm)}")
if len(flap_comm) > 0:
    print("Most efficient flap communication:")
    print(flap_comm.head())

# Find flap-to-active site communication
flap_active_comm = high_eff[
    ((high_eff['Source_Residue'].isin(flap_a + flap_b)) & (high_eff['Target_Residue'].isin(active_site))) |
    ((high_eff['Source_Residue'].isin(active_site)) & (high_eff['Target_Residue'].isin(flap_a + flap_b)))
]

print(f"Flap-to-active site communication: {len(flap_active_comm)}")
```

### **2. Network Properties Validation**
```python
# Load network metrics
import json
with open('test_fixed_analysis_results/01_residue_interaction_networks/test_fixed_analysis_community_analysis.json', 'r') as f:
    community_data = json.load(f)

print("Network Community Analysis:")
print(f"  Method: {community_data.get('method', 'N/A')}")
print(f"  Communities: {community_data.get('n_communities', 'N/A')}")
print(f"  Modularity: {community_data.get('modularity', 'N/A'):.4f}")

# Check centrality significance
centrality_sig = pd.read_csv('test_fixed_analysis_results/01_residue_interaction_networks/test_fixed_analysis_significant_nodes.csv')
print(f"Statistically significant nodes: {len(centrality_sig)}")
```

## **Performance Validation**

### **1. Analysis Timing**
```bash
# Time the analysis
time python md_compare_cli.py single \
  -t your_protein.pdb \
  -x your_trajectory.xtc \
  -n performance_test \
  --compute-msm
```

### **2. Memory Usage Monitoring**
```python
# Monitor memory during analysis
import psutil
import subprocess
import time

def monitor_analysis():
    process = subprocess.Popen([
        'python', 'md_compare_cli.py', 'single',
        '-t', 'your_protein.pdb',
        '-x', 'your_trajectory.xtc', 
        '-n', 'memory_test'
    ])
    
    max_memory = 0
    while process.poll() is None:
        try:
            memory_mb = psutil.Process(process.pid).memory_info().rss / 1024 / 1024
            max_memory = max(max_memory, memory_mb)
            time.sleep(1)
        except psutil.NoSuchProcess:
            break
    
    print(f"Maximum memory usage: {max_memory:.1f} MB")
    return process.returncode == 0

# Run memory monitoring
success = monitor_analysis()
print(f"Analysis completed successfully: {success}")
```

## **Success Criteria**

Your implementation is working correctly if:

### ✅ **Functional Requirements**
- [ ] Analysis completes without errors
- [ ] All 7 subdirectories are created  
- [ ] Communication efficiency matrix is 198×198 (all residue pairs)
- [ ] High-efficiency pairs CSV contains meaningful data
- [ ] MSM transition matrix is properly normalized (row sums ≈ 1.0)
- [ ] Energy landscape visualization is generated (simple or detailed)

### ✅ **Data Quality Requirements**  
- [ ] Communication matrix has non-zero diagonal (perfect self-communication)
- [ ] High-efficiency pairs show cross-chain communication for HIV protease
- [ ] MSM shows reasonable connectivity (>80% connected states)
- [ ] Implied timescales show hierarchical separation
- [ ] Network communities are detected with reasonable modularity (>0.3)

### ✅ **Output Quality Requirements**
- [ ] All CSV files are properly formatted with headers
- [ ] Visualizations are high-resolution (300 DPI)
- [ ] Summary reports contain meaningful scientific information
- [ ] File organization is clear and professional

## 🔧 **Troubleshooting**

If any tests fail, check:

1. **File Paths**: Ensure input files exist and paths are correct
2. **Dependencies**: Verify PyEMMA, pandas, matplotlib are installed
3. **Memory**: Large systems may need more RAM or stride parameters
4. **Permissions**: Ensure write access to output directory
5. **Trajectory Format**: Verify MDAnalysis can read your trajectory format

## **Excel Analysis Guide**

Open the communication efficiency matrix in Excel:

1. **Load Matrix**: Open `*_communication_efficiency_matrix.csv`
2. **Create Heatmap**: Select data → Insert → Charts → Heatmap
3. **Filter High Values**: Use conditional formatting to highlight efficiency >0.3
4. **Cross-Chain Analysis**: Filter rows/columns by chain (A_ vs B_ residues)
5. **Pathway Identification**: Use `*_high_efficiency_pairs.csv` for pathway analysis

---

**Testing Complete!** ✅ Your enhanced MD-Compare v1.4.0 should now provide comprehensive protein dynamics analysis with full communication efficiency matrices and organized output structure! 🧬📊✨
