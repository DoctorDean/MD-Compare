def build_markov_state_model(universe, sub_traj_ranges, cv_selections):
    """Build Markov State Models for kinetic analysis"""
    msm_results = {}

    for i, (start, end) in enumerate(sub_traj_ranges):
        # Calculate features (collective variables)
        features = []
        for cv_def in cv_selections:
            cv_values = calculate_collective_variable(universe, cv_def, start, end)
            features.append(cv_values)

        features = np.array(features).T

        # Cluster conformational space
        clustering = pyemma.coordinates.cluster_kmeans(features, k=100, max_iter=500)
        discretized_trajectory = clustering.dtrajs[0]

        # Build MSM
        lag_time = 10  # frames
        msm = pyemma.msm.estimate_markov_model(discretized_trajectory, lag=lag_time)

        msm_results[f'sub_traj_{i}'] = {
            'msm': msm,
            'clustering': clustering,
            'stationary_probabilities': msm.stationary_distribution,
            'transition_matrix': msm.transition_matrix,
            'implied_timescales': msm.timescales()
        }

    return msm_results