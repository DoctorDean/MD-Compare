from MDContactNetworks import ContactNetwork
import networkx as nx


def analyze_residue_networks(universe, sub_traj_ranges):
    """Analyze residue interaction networks"""
    network_results = {}

    for i, (start, end) in enumerate(sub_traj_ranges):
        # Create contact network
        cn = ContactNetwork(universe,
                            selection1="protein",
                            selection2="protein",
                            cutoff=4.5,  # Å
                            start=start,
                            stop=end)

        # Calculate network properties
        G = cn.graph

        network_metrics = {
            'betweenness_centrality': nx.betweenness_centrality(G),
            'closeness_centrality': nx.closeness_centrality(G),
            'clustering_coefficient': nx.clustering(G),
            'shortest_path_lengths': dict(nx.all_pairs_shortest_path_length(G))
        }

        # Community detection
        communities = nx.community.greedy_modularity_communities(G)

        network_results[f'sub_traj_{i}'] = {
            'metrics': network_metrics,
            'communities': communities,
            'persistent_contacts': cn.get_persistent_contacts(threshold=0.7)
        }

    return network_results


def map_allosteric_pathways(network_results, mutation_sites, active_site):
    """Map allosteric communication pathways"""
    pathways = {}

    for sub_traj, data in network_results.items():
        G = data['graph']
        paths = {}

        for mut_site in mutation_sites:
            for active_res in active_site:
                try:
                    path = nx.shortest_path(G, mut_site, active_res)
                    paths[f'{mut_site}_to_{active_res}'] = path
                except nx.NetworkXNoPath:
                    paths[f'{mut_site}_to_{active_res}'] = None

        pathways[sub_traj] = paths

    return pathways