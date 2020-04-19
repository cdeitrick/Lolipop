import itertools
from typing import *

import numpy
import pandas
from loguru import logger
from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller.clustering.metrics.distance_cache import DistanceCache
	from muller.dataio import projectdata
except ModuleNotFoundError:
	from muller.clustering.metrics.distance_cache import DistanceCache
	from muller.dataio import projectdata


def format_linkage_matrix(linkage_table, total_members: Optional[int]) -> pandas.DataFrame:
	linkage_dataframe = pandas.DataFrame(linkage_table, columns = ["left", "right", "distance", "observations"])

	linkage_dataframe['left'] = linkage_dataframe['left'].astype(int)
	linkage_dataframe['right'] = linkage_dataframe['right'].astype(int)
	linkage_dataframe['observations'] = linkage_dataframe['observations'].astype(int)
	if total_members:
		linkage_dataframe['clusterId'] = list(total_members + i for i in range(len(linkage_dataframe)))
	# Rearrange the columns to be more readable.
	linkage_dataframe = linkage_dataframe[['left', 'right', 'distance', 'observations', 'clusterId']]
	return linkage_dataframe


class HierarchalCluster:
	def __init__(self, linkage: str = 'ward', cluster: str = 'distance'):
		self.linkage_method = linkage
		self.cluster_method = cluster

	@staticmethod
	def _add_starting_genotypes(pair_array, starting_genotypes) -> Dict[str, str]:
		for genotype in starting_genotypes:
			combinations = itertools.combinations(genotype, 2)
			values = dict()
			for a, b in combinations:
				values[a, b] = 0
				values[b, a] = 0
			pair_array = pair_array.update(values)
		return pair_array

	@staticmethod
	def _cluster_by_distance(linkage_table: numpy.ndarray, cutoff: float) -> numpy.ndarray:
		""" Try to infer a good distance cutoff by detecting the first changepoint in the sorted array of distances."""
		clusters = hierarchy.fcluster(linkage_table, t = cutoff, criterion = 'distance')
		return clusters

	@staticmethod
	def _cluster_by_monocrit(linkage_table: numpy.ndarray, cutoff: float, inconsistent: pandas.DataFrame) -> numpy.ndarray:
		MR = hierarchy.maxRstat(linkage_table, inconsistent.values, 1)
		clusters = hierarchy.fcluster(linkage_table, t = cutoff, criterion = 'monocrit', monocrit = MR)
		return clusters

	@staticmethod
	def _cluster_by_inconsistent(linkage_table: numpy.ndarray, cutoff: float, inconsistent: pandas.DataFrame) -> numpy.ndarray:
		clusters = hierarchy.fcluster(linkage_table, t = cutoff, criterion = 'inconsistent', R = inconsistent.values)
		return clusters

	@staticmethod
	def _get_inconsistent(distances: numpy.ndarray) -> pandas.DataFrame:
		inconsistent = hierarchy.inconsistent(distances, 10)
		inc = pandas.DataFrame(inconsistent)
		inc.columns = ["mean", "std", "count", "coefficient"]
		return inc

	@staticmethod
	def _label_clusters(clusters: numpy.ndarray, labels: Iterable[str]) -> List[List[str]]:
		""" Maps the original labels to the corresponding numerical index in `clusters`"""
		cluster_map = {i: list() for i in clusters}
		for i, j in zip(clusters, labels):
			cluster_map[i].append(j)
		return list(cluster_map.values())


	def link_clusters(self, distances: numpy.ndarray, num: int) -> pandas.DataFrame:
		Z = hierarchy.linkage(distances, method = self.linkage_method, optimal_ordering = True)

		return format_linkage_matrix(Z, num)

	def cluster(self, linkage_table: pandas.DataFrame, cutoff: float, labels: Optional[Iterable[str]] = None) -> List[List[Union[int, str]]]:
		inconsistent = self._get_inconsistent(linkage_table.values)
		if self.cluster_method == 'distance':
			clusters = self._cluster_by_distance(linkage_table.values, cutoff)
		elif self.cluster_method == 'monocrit':
			clusters = self._cluster_by_monocrit(linkage_table.values, cutoff, inconsistent)

		elif self.cluster_method == 'inconsistent':
			clusters = self._cluster_by_inconsistent(linkage_table.values, cutoff, inconsistent)
		else:
			message = f"'{self.cluster_method}' is not a currently supported clustering method when using hierarchical clustering."
			logger.error(message)
			raise ValueError(message)

		if labels is not None:
			clusters = self._label_clusters(clusters, labels)

		return clusters
	@staticmethod
	def adjust_similarity_cutoff(quantile:float, distances: List[float])->float:
		""" Adjusts the `similarity_cutoff` value to work with the distance observations."""

		#result = pandas.Series(distances).quantile(quantile)
		result = pandas.Series([i for i in distances if 0 < i < max(distances)]).quantile(quantile)
		return result


	def run(self, pair_array: DistanceCache, starting_genotypes: List[List[str]] = None, similarity_cutoff: Optional[float] = None) -> projectdata.DataHierarchalCluster:
		"""
		Parameters
		----------
		pair_array
		starting_genotypes: List[List[str]]
			Each element should be a list of trajectories known to be in the same genotype.
		similarity_cutoff: Optional[float]
			If not given, the similarity cutoff will be generated automatically.
		"""

		# If known genotypes are given, modify the pair_array so that they will be grouped together.
		if starting_genotypes:
			pair_array = self._add_starting_genotypes(pair_array, starting_genotypes)
		squaremap = pair_array.squareform()
		distance_array = distance.squareform(squaremap.values)
		linkage_table = self.link_clusters(distance_array, len(squaremap.index))
		reduced_linkage_table = linkage_table[['left', 'right', 'distance', 'observations']]  # Removes the extra column

		if similarity_cutoff is None:
			quantile = 0.05
		else:
			quantile = similarity_cutoff

		distance_cutoff = self.adjust_similarity_cutoff(quantile, list(pair_array.pairwise_values.values()))

		logger.debug(f"Using Hierarchical Clustering with similarity cutoff {distance_cutoff}")

		clusters = self.cluster(reduced_linkage_table, distance_cutoff, squaremap.index)

		result = projectdata.DataHierarchalCluster(
			clusters = clusters,
			table_linkage = linkage_table,
			distance_cutoff = distance_cutoff,
			distance_quantile = quantile
		)

		return result

