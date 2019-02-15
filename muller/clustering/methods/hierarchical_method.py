from typing import Any, List, Tuple
import pandas
from scipy.cluster import hierarchy
from scipy.spatial import distance
import logging
logger = logging.getLogger(__file__)
try:
	from clustering.metrics.pairwise_calculation_cache import PairwiseCalculationCache
except ModuleNotFoundError:
	from ..metrics.pairwise_calculation_cache import PairwiseCalculationCache


def hierarchical_method(pair_array: PairwiseCalculationCache, similarity_cutoff: float, cluster_method: str = 'monocrit') -> Tuple[List[List[str]], Any]:
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
	logger.debug(f"similarity_cutoff: {similarity_cutoff}")
	logger.debug(f"cluster_method: {cluster_method}")
	squaremap = pair_array.squareform()
	condensed_squaremap = distance.squareform(squaremap.values)

	Z = hierarchy.linkage(condensed_squaremap, method = 'ward', optimal_ordering = True)
	inconsistent = hierarchy.inconsistent(Z, 10)

	MR = hierarchy.maxRstat(Z, inconsistent, 1)

	if cluster_method == 'distance':
		maximum_distance = max(pair_array.pairwise_values.values())
		distances = pandas.Series([i for i in pair_array.pairwise_values.values() if (i != maximum_distance and i > 0)])
		similarity_cutoff = distances.quantile(.1)

		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'distance')
	elif cluster_method == 'monocrit':
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'monocrit', monocrit = MR)
	elif cluster_method == 'inconsistent':
		clusters = hierarchy.fcluster(Z, t = similarity_cutoff, criterion = 'inconsistent')
	else:
		message = f"'{cluster_method}' is not a currently supported clustering method when using hierarchical clustering."
		logger.error(message)
		raise ValueError(message)
	cluster_map = {i: list() for i in clusters}
	for i, j in zip(clusters, squaremap.index):
		cluster_map[i].append(j)
	return list(cluster_map.values()), Z


if __name__ == "__main__":
	pass
