#!/usr/bin/env python3
"""
MD-Compare: Comprehensive Protein Dynamics Analysis Platform

A comprehensive toolkit for analyzing molecular dynamics simulations with advanced 
network analysis, conformational dynamics, and kinetic modeling capabilities.

Version: 1.4.0
Author: Dr Dean Sherry
License: MIT

Key Features:
- Network topology analysis with community detection
- Dynamic cross-correlation matrix (DCCM) analysis  
- Principal component analysis (PCA) and energy landscapes
- Markov State Model (MSM) analysis with PyEMMA integration
- Allosteric pathway mapping and hotspot identification
- Statistical validation and cross-platform compatibility

Example Usage:
    >>> from md_compare_core import NetworkAnalyzer, AnalysisConfig
    >>> config = AnalysisConfig(compute_msm=True, msm_lag_time=10)
    >>> analyzer = NetworkAnalyzer(config)
    >>> # ... perform analysis
"""

__version__ = "1.4.0"
__author__ = "Dr Dean Sherry"
__email__ = ""
__license__ = "MIT"
__url__ = "https://github.com/DoctorDean/md-compare"

# Version info tuple for programmatic access
__version_info__ = tuple(int(v) for v in __version__.split('.'))

# Import main classes for convenient access
try:
    from .md_compare_core import (
        NetworkAnalyzer,
        AnalysisConfig,
        NetworkMetrics,
        MDSimulation,
        OutputManager
    )
    
    # Make key classes available at package level
    __all__ = [
        'NetworkAnalyzer',
        'AnalysisConfig', 
        'NetworkMetrics',
        'MDSimulation',
        'OutputManager',
        '__version__',
        '__version_info__'
    ]
    
except ImportError:
    # Fallback for development/testing scenarios
    __all__ = ['__version__', '__version_info__']

# Optional dependency checks and feature flags
FEATURES = {
    'pyemma': False,
    'igraph': False,
    'leidenalg': False,
    'pandas': False,
}

try:
    import pyemma
    FEATURES['pyemma'] = True
except ImportError:
    pass

try:
    import igraph
    FEATURES['igraph'] = True
except ImportError:
    pass

try:
    import leidenalg
    FEATURES['leidenalg'] = True
except ImportError:
    pass

try:
    import pandas
    FEATURES['pandas'] = True
except ImportError:
    pass

def check_dependencies(verbose=False):
    """
    Check availability of optional dependencies
    
    Parameters:
    -----------
    verbose : bool
        If True, print detailed dependency information
        
    Returns:
    --------
    dict : Feature availability status
    """
    if verbose:
        print(f"MD-Compare v{__version__}")
        print("=" * 40)
        print("Core Dependencies:")
        
        # Check core dependencies
        core_deps = ['MDAnalysis', 'networkx', 'numpy', 'scipy', 
                    'matplotlib', 'seaborn', 'scikit-learn']
        
        for dep in core_deps:
            try:
                __import__(dep.lower().replace('-', '_'))
                print(f"  ✓ {dep}")
            except ImportError:
                print(f"  ✗ {dep} (REQUIRED)")
        
        print("\nOptional Dependencies:")
        for feature, available in FEATURES.items():
            status = "✓" if available else "✗"
            print(f"  {status} {feature}")
            
        if FEATURES['pyemma']:
            print("    → Markov State Model analysis available")
        if FEATURES['igraph'] and FEATURES['leidenalg']:
            print("    → Advanced community detection available")
        if FEATURES['pandas']:
            print("    → Excel export functionality available")
    
    return FEATURES

def get_version_info():
    """Get comprehensive version and dependency information"""
    import sys
    
    info = {
        'md_compare_version': __version__,
        'python_version': sys.version,
        'platform': sys.platform,
        'features': FEATURES,
    }
    
    # Add version info for key dependencies
    try:
        import MDAnalysis
        info['mdanalysis_version'] = MDAnalysis.__version__
    except (ImportError, AttributeError):
        info['mdanalysis_version'] = 'Not available'
    
    try:
        import networkx
        info['networkx_version'] = networkx.__version__
    except (ImportError, AttributeError):
        info['networkx_version'] = 'Not available'
    
    if FEATURES['pyemma']:
        try:
            import pyemma
            info['pyemma_version'] = pyemma.__version__
        except AttributeError:
            info['pyemma_version'] = 'Unknown'
    
    return info

# Package metadata
__doc__ += f"""
Version Information:
    MD-Compare: {__version__}
    
Available Features:
    PyEMMA (MSM analysis): {FEATURES['pyemma']}
    Advanced communities: {FEATURES['igraph'] and FEATURES['leidenalg']}
    Excel export: {FEATURES['pandas']}

For dependency checking, use:
    >>> import md_compare
    >>> md_compare.check_dependencies(verbose=True)
"""
