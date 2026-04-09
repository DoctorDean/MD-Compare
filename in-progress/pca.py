from sklearn.decomposition import PCA
from MDAnalysis.analysis import pca


def perform_pca_analysis(universe, sub_traj_ranges):
    """Perform PCA on conformational dynamics"""
    pca_results = {}

    for i, (start, end) in enumerate(sub_traj_ranges):
        # PCA analysis
        pc_analysis = pca.PCA(universe,
                              select="name CA",
                              align=True,
                              start=start,
                              stop=end)
        pc_analysis.run()

        # Extract results
        eigenvalues = pc_analysis.eigenvalues
        eigenvectors = pc_analysis.p_components

        # Project trajectory onto PC space
        transformed = pc_analysis.transform(universe.select_atoms("name CA"),
                                            start=start, stop=end)

        pca_results[f'sub_traj_{i}'] = {
            'eigenvalues': eigenvalues,
            'eigenvectors': eigenvectors,
            'transformed_trajectory': transformed,
            'cumulative_variance': np.cumsum(eigenvalues) / np.sum(eigenvalues)
        }

    return pca_results


def difference_pca(wt_pca_results, mutant_pca_results):
    """Perform difference PCA between WT and mutant"""
    # Compare essential subspaces
    subspace_overlap = {}

    for i in range(5):  # For each sub-trajectory
        wt_evecs = wt_pca_results[f'sub_traj_{i}']['eigenvectors'][:, :10]  # First 10 PCs
        mut_evecs = mutant_pca_results[f'sub_traj_{i}']['eigenvectors'][:, :10]

        # Calculate subspace overlap
        overlap_matrix = np.dot(wt_evecs.T, mut_evecs)
        subspace_overlap[f'sub_traj_{i}'] = np.sum(overlap_matrix ** 2)

    return subspace_overlap