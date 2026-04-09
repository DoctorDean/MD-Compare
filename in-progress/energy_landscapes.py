import pyemma


def construct_free_energy_landscapes(universe, sub_traj_ranges, cv_selections):
    """Construct 2D free energy landscapes"""
    landscape_results = {}

    for i, (start, end) in enumerate(sub_traj_ranges):
        # Calculate collective variables
        cvs = []
        for cv_def in cv_selections:
            cv_values = calculate_collective_variable(universe, cv_def, start, end)
            cvs.append(cv_values)

        cvs = np.array(cvs).T  # Shape: (n_frames, n_cvs)

        # Estimate free energy using histograms
        hist, edges = np.histogramdd(cvs, bins=50)
        hist[hist == 0] = 1  # Avoid log(0)

        # Convert to free energy (kT units)
        free_energy = -np.log(hist)
        free_energy -= np.min(free_energy)  # Set minimum to 0

        landscape_results[f'sub_traj_{i}'] = {
            'free_energy': free_energy,
            'edges': edges,
            'cvs': cvs
        }

    return landscape_results


def calculate_collective_variable(universe, cv_definition, start, end):
    """Calculate collective variables (e.g., flap distance, pocket volume)"""
    cv_values = []

    for ts in universe.trajectory[start:end]:
        if cv_definition['type'] == 'distance':
            # Example: flap tip distance
            atom1 = universe.select_atoms(cv_definition['selection1'])
            atom2 = universe.select_atoms(cv_definition['selection2'])
            distance = np.linalg.norm(atom1.center_of_mass() - atom2.center_of_mass())
            cv_values.append(distance)
        elif cv_definition['type'] == 'volume':
            # Pocket volume calculation
            volume = calculate_pocket_volume_frame(universe, cv_definition['pocket_residues'])
            cv_values.append(volume)

    return np.array(cv_values)