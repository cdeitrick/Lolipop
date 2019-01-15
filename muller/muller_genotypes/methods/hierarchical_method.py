from typing import Any, List, Tuple

from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller_genotypes.metrics.pairwise_calculation import PairwiseCalculation
except ModuleNotFoundError:
	from ..metrics.pairwise_calculation import PairwiseCalculation

def hierarchical_method(pair_array: PairwiseCalculation, similarity_cutoff: float, cluster_method:str = 'monocrit') -> Tuple[List[List[str]], Any]:
	"""

	Parameters
	----------
	pair_array
	similarity_cutoff
	cluster_method: {'distance;, 'monocrit'}; default 'monocrit'
		The method to cluster leafs by.
	Returns
	-------

	"""
	squaremap = pair_array.squareform('X')
	condensed_squaremap = distance.squareform(squaremap.values)

	Z = hierarchy.linkage(condensed_squaremap)

	import matplotlib.pyplot as plt
	plt.figure(figsize = (10, 10))
	plt.title('Hierarchical Clustering Dendrogram')
	plt.xlabel('sample index')
	plt.ylabel('distance')

	hierarchy.dendrogram(
		Z,
		color_threshold = similarity_cutoff,
		leaf_rotation = 90,  # rotates the x axis labels
		leaf_font_size = 8,  # font size for the x axis labels,
		labels = squaremap.index
	)
	plt.savefig("clusterimage.png")
	if cluster_method == 'distance':
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'distance')
	elif cluster_method == 'monocrit':
		inconsistent = hierarchy.inconsistent(Z)
		MR = hierarchy.maxRstat(Z, inconsistent, 1)
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'monocrit', monocrit = MR)
	else:
		message = f"'{cluster_method}' is not a currently supported clustering method when using hierarchical clustering."
		raise ValueError(message)
	cluster_map = {i: list() for i in clusters}
	for i, j in zip(clusters, squaremap.index):
		cluster_map[i].append(j)
	return list(cluster_map.values()), Z
