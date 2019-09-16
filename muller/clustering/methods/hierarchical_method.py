import itertools
from typing import Any, List, Tuple

import pandas
from loguru import logger
from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller.clustering.metrics.distance_cache import DistanceCache
except ModuleNotFoundError:
	from ..metrics.distance_cache import DistanceCache


class AgglomerativeClustering:
	"""
	Parameters
	----------
	pair_array
	similarity_cutoff
	cluster_method: {'distance;, 'monocrit'}; default 'monocrit'
		The method to cluster leafs by.
	starting_genotypes: List[List[str]]
		Each element should be a list of trajectories known to be in the same genotype.
	"""

	def __init__(self, similarity_cutoff: float, cluster_method: str = 'distance'):
		self.similarity_cutoff = similarity_cutoff,
		self.cluster_method = cluster_method
	@staticmethod
	def apply_starting_genotypes(pair_array: DistanceCache, starting_genotypes: List[List[str]] = None) -> DistanceCache:
		""" Modified the distances in the pair array so that each of the trajectories grouped as genotypes are forced together."""
		for genotype in starting_genotypes:
			combinations = itertools.combinations(genotype, 2)
			values = dict()
			for a, b in combinations:
				values[a, b] = 0
				values[b, a] = 0
			pair_array = pair_array.update(values)
		return pair_array

	def run(self, pair_array: DistanceCache, starting_genotypes: List[List[str]] = None, linkage_method:str = 'complete', cluster_method: str = 'distance'):
		if starting_genotypes:
			pair_array = self.apply_starting_genotypes(pair_array, starting_genotypes)

		squaremap = pair_array.squareform()
		condensed_squaremap = distance.squareform(squaremap.values)
		Z = hierarchy.linkage(condensed_squaremap, method = linkage_method, optimal_ordering = True)
		inconsistent = hierarchy.inconsistent(Z, 10)
		MR = hierarchy.maxRstat(Z, inconsistent, 1)

	def between_cluster_variation(self, distances:DistanceCache, genotypes:List[List[str]]):
		import itertools
		avarage = sum(distances.values()) / len(distances)

		genotype_distances = dict()
		for genotype in genotypes:
			number_of_points = len(genotype)
			all_possible_combinations = itertools.combinations(genotype, 2)
			total = 0
			for i,j in all_possible_combinations:
				distance = distances.get(i,j)
				total += distance

	def get_within_cluster_variation(self, cluster:List[str], distances:DistanceCache)->float:
		combs = itertools.combinations(cluster, 2)
		total = sum(distances.get(i,j) for i,j in combs)

	def get_cluster_variation(self, cluster:List[str], distances:DistanceCache):
		combs = itertools.combinations(cluster, 2)
		total = sum(distances.get(i,j) for i,j in combs)



def hierarchical_method(pair_array: DistanceCache, similarity_cutoff: float, cluster_method: str = 'distance',
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
	z = pandas.DataFrame(Z)
	z.columns = ['left', 'right', 'distance', 'count']
	z.to_csv("z.tsv", sep = "\t", index = False)
	inconsistent = hierarchy.inconsistent(Z, 10)
	MR = hierarchy.maxRstat(Z, inconsistent, 1)
	inc = pandas.DataFrame(inconsistent)
	inc.columns = ["mean", "std", "count", "coefficient"]
	# inc.to_csv(f"inconsistent.tsv", sep = "\t", index = False)
	logger.info(inc['std'].values)
	similarity_cutoff = inc['std'].iloc[1] * 4

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
