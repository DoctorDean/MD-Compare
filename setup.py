#!/usr/bin/env python3

"""
Setup script for MD-Compare: Molecular Dynamics Simulation Comparison Toolkit
"""

import sys
from pathlib import Path
from setuptools import setup, find_packages

# Ensure Python version compatibility
if sys.version_info < (3, 8):
    raise RuntimeError("This package requires Python 3.8 or later")

# Read README for long description
readme_path = Path(__file__).parent / "README.md"
if readme_path.exists():
    with open(readme_path, "r", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "MD-Compare: Molecular Dynamics Simulation Comparison Toolkit"

# Read requirements
requirements_path = Path(__file__).parent / "requirements.txt"
if requirements_path.exists():
    with open(requirements_path, "r", encoding="utf-8") as f:
        requirements = [
            line.strip() 
            for line in f 
            if line.strip() and not line.startswith("#")
        ]
else:
    requirements = [
        "MDAnalysis>=2.4.0",
        "networkx>=2.8",
        "numpy>=1.21.0", 
        "scipy>=1.7.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "pandas>=1.3.0"
    ]

# Optional dependencies
extras_require = {
    'dev': [
        'pytest>=6.0.0',
        'pytest-cov>=3.0.0', 
        'black>=22.0.0',
        'flake8>=4.0.0',
        'sphinx>=4.0.0',
        'sphinx-rtd-theme>=1.0.0'
    ],
    'viz': [
        'plotly>=5.0.0',
        'ipywidgets>=7.6.0',
        'pygraphviz>=1.7'
    ],
    'ml': [
        'scikit-learn>=1.0.0'
    ],
    'notebook': [
        'jupyter>=1.0.0',
        'ipykernel>=6.0.0'
    ]
}

# All extras combined
extras_require['all'] = [
    dep for deps in extras_require.values() for dep in deps
]

setup(
    name="md-compare",
    version="1.0.0",
    description="Comprehensive toolkit for comparing molecular dynamics simulations using network analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    
    author="Your Name", 
    author_email="your.email@institution.edu",
    url="https://github.com/yourusername/md-compare",
    
    packages=find_packages(),
    
    # Main scripts
    scripts=[
        "md_compare_cli.py",
        "md_compare_core.py",
        "utils.py"
    ],
    
    # Console entry points
    entry_points={
        'console_scripts': [
            'md-compare=md_compare_cli:main',
            'md-compare-core=md_compare_core:main',
        ],
    },
    
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require=extras_require,
    
    # Package metadata
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research", 
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9", 
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    
    keywords="molecular dynamics, network analysis, protein dynamics, simulation comparison, bioinformatics, computational biology",
    
    # Include additional files
    include_package_data=True,
    zip_safe=False,
    
    # Project URLs
    project_urls={
        "Bug Reports": "https://github.com/yourusername/md-compare/issues",
        "Source": "https://github.com/yourusername/md-compare",
        "Documentation": "https://github.com/yourusername/md-compare/wiki",
    },
)

