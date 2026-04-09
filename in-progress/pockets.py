def analyze_binding_pocket_dynamics(universe, sub_traj_ranges, pocket_residues):
    """Comprehensive binding pocket analysis"""
    pocket_results = {}

    for i, (start, end) in enumerate(sub_traj_ranges):
        pocket_analysis = {
            'volumes': [],
            'surface_areas': [],
            'sphericity': [],
            'hydrophobicity': [],
            'electrostatic_potential': []
        }

        for ts in universe.trajectory[start:end]:
            # Pocket volume using Alpha shapes
            volume = calculate_alpha_shape_volume(universe, pocket_residues)
            pocket_analysis['volumes'].append(volume)

            # Surface area
            surface_area = calculate_pocket_surface_area(universe, pocket_residues)
            pocket_analysis['surface_areas'].append(surface_area)

            # Shape descriptors
            sphericity = calculate_pocket_sphericity(universe, pocket_residues)
            pocket_analysis['sphericity'].append(sphericity)

        pocket_results[f'sub_traj_{i}'] = pocket_analysis

    return pocket_results


def calculate_alpha_shape_volume(universe, pocket_residues):
    """Calculate pocket volume using alpha shapes"""
    from scipy.spatial import ConvexHull
    import alphashape

    # Get coordinates of pocket-lining atoms
    pocket_atoms = universe.select_atoms(f"resid {' '.join(map(str, pocket_residues))} and not name H*")
    coords = pocket_atoms.positions

    # Calculate alpha shape
    alpha_shape = alphashape.alphashape(coords, 2.0)

    # Estimate volume (simplified - use convex hull as approximation)
    try:
        hull = ConvexHull(coords)
        return hull.volume
    except:
        return 0.0


def run_fpocket_analysis(pdb_file, output_dir):
    """Run fpocket analysis on representative structures"""
    import subprocess

    cmd = f"fpocket -f {pdb_file} -d {output_dir}"
    subprocess.run(cmd, shell=True)

    # Parse fpocket output
    pocket_file = f"{output_dir}/pockets/pocket1_atm.pdb"
    # Extract pocket properties from fpocket output files

    return parse_fpocket_results(output_dir)