from typing import Any, List, Tuple

from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller_genotypes.metrics.pairwise_calculation_cache import PairwiseCalculation
except ModuleNotFoundError:
	from ..metrics.pairwise_calculation_cache import PairwiseCalculation


def hierarchical_method(pair_array: PairwiseCalculation, similarity_cutoff: float, cluster_method: str = 'monocrit') -> Tuple[List[List[str]], Any]:
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
	squaremap = pair_array.squareform()
	condensed_squaremap = distance.squareform(squaremap.values)

	Z = hierarchy.linkage(condensed_squaremap, method = 'ward')
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


if __name__ == "__main__":
	pass