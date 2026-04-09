def calculate_cross_correlations(universe, sub_traj_ranges):
    """Calculate dynamic cross-correlation matrices"""
    from MDAnalysis.analysis import align

    correlation_results = {}

    for i, (start, end) in enumerate(sub_traj_ranges):
        # Align trajectory
        aligner = align.AlignTraj(universe, universe,
                                  select="name CA",
                                  start=start, stop=end)
        aligner.run()

        # Calculate position deviations
        ca_atoms = universe.select_atoms("name CA")
        n_residues = len(ca_atoms)
        n_frames = end - start

        positions = np.zeros((n_frames, n_residues, 3))

        for j, ts in enumerate(universe.trajectory[start:end]):
            positions[j] = ca_atoms.positions

        # Calculate correlations
        mean_pos = np.mean(positions, axis=0)
        deviations = positions - mean_pos

        correlation_matrix = np.zeros((n_residues, n_residues))

        for res1 in range(n_residues):
            for res2 in range(n_residues):
                dev1 = deviations[:, res1, :].flatten()
                dev2 = deviations[:, res2, :].flatten()

                correlation_matrix[res1, res2] = np.corrcoef(dev1, dev2)[0, 1]

        correlation_results[f'sub_traj_{i}'] = correlation_matrix

    return correlation_results









#Draw the graph
file=np.loadtxt(working_dir + f"cross_corr_{system_name}.csv")
data= np.array(file)

fig, ax = plt.subplots(figsize=(11,9))

# Set the vmin and vmax values according to your interests
# You can change the cmap styles (PiYG, PRGn, BrBG, PuOr, RdGy, RdBu,
# RdYlBu,RdYlGn, Spectral, coolwarm, bwr, seismic, twilight)
im = ax.imshow(data, cmap=plt.cm.viridis, vmin=-1, vmax=1, origin='lower')
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.tick_params(labelsize=24)
for tick in cbar.ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
ax.set_title(f"Cross Correlation {system_name}", fontname= "Times New Roman", fontsize=42, pad=15)
ax.tick_params(labelsize=24)
for tick in ax.get_xticklabels():
    tick.set_fontname("Times New Roman")
for tick in ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
fig.tight_layout()


# Set major tick locations (e.g., every 2 units)
ax.set_xticks(np.arange(0, data.shape[1], 100))
ax.set_yticks(np.arange(0, data.shape[0], 100))

# Set minor tick locations (e.g., every 1 unit)
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(MultipleLocator(50))

# Enable minor ticks
ax.tick_params(which='minor', direction='out', length=4, color='gray')

#Save the graph
fig.savefig(working_dir + f'cross_corr_{system_name}.jpeg', dpi=500)