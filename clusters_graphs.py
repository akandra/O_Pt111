import numpy as np
from scipy.spatial import distance_matrix
import networkx as nx

# Example: your lattice points as Nx2 array
points = np.array([
    [0, 0],
    [1, 0],
    [0.5, np.sqrt(3)/2],
    [2, 0],
    [3, 0],
    [1.5, np.sqrt(3)/2],
    # add more points
])

# Define your third-nearest-neighbor distance
# On a hexagonal lattice with unit spacing:
# 1st NN ≈ 1.0
# 2nd NN ≈ sqrt(3)
# 3rd NN ≈ 2.0
cutoff = 2.0 + 1e-6

# Compute pairwise distances
dist = distance_matrix(points, points)

# Create graph: connect points within cutoff
G = nx.Graph()
for i in range(len(points)):
    for j in range(i + 1, len(points)):
        if dist[i, j] <= cutoff:
            G.add_edge(i, j)

# Find connected components (clusters)
clusters = list(nx.connected_components(G))

# Print results
for idx, cluster in enumerate(clusters):
    print(f"Cluster {idx + 1}: {sorted(cluster)}")

    # this is a very unimportant change.
