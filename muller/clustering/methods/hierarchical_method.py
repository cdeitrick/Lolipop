import itertools
from typing import Any, List, Tuple

import pandas
from loguru import logger
from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller.clustering.metrics.pairwise_calculation_cache import PairwiseCalculationCache
except ModuleNotFoundError:
	from ..metrics.pairwise_calculation_cache import PairwiseCalculationCache


def hierarchical_method(pair_array: PairwiseCalculationCache, similarity_cutoff: float, cluster_method: str = 'distance',
		starting_genotypes: List[List[str]] = None) -> Tuple[List[List[str]], Any]:
	"""

	Parameters
	----------
	pair_array
	similarity_cutoff
	cluster_method: {'distance;, 'monocrit'}; default 'monocrit'
		The method to cluster leafs by.
	starting_genotypes: List[List[str]]
		Each element should be a list of trajectories known to be in the same genotype.
	Returns
	-------

	"""
	logger.debug(f"similarity_cutoff: {similarity_cutoff}")
	logger.debug(f"cluster_method: {cluster_method}")

	# If known genotypes are given, modify the pair_array so that they will be grouped together.
	if starting_genotypes:
		for genotype in starting_genotypes:
			combinations = itertools.combinations(genotype, 2)
			values = dict()
			for a, b in combinations:
				values[a, b] = 0
				values[b, a] = 0
			pair_array = pair_array.update(values)

	squaremap = pair_array.squareform()
	condensed_squaremap = distance.squareform(squaremap.values)
	method = 'complete'
	Z = hierarchy.linkage(condensed_squaremap, method = method, optimal_ordering = True)
	inconsistent = hierarchy.inconsistent(Z, 10)
	MR = hierarchy.maxRstat(Z, inconsistent, 1)

	if cluster_method == 'distance':
		maximum_distance = max(pair_array.pairwise_values.values())
		distances = pandas.Series([i for i in pair_array.pairwise_values.values() if (i != maximum_distance and i > 0)])

		similarity_cutoff = distances.quantile(similarity_cutoff)

		logger.debug(f"Using Hierarchical Clustering with similarity cutoff {similarity_cutoff}")

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
