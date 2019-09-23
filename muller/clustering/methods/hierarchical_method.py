import itertools
import statistics
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy
import pandas
from loguru import logger
from scipy.cluster import hierarchy
from scipy.spatial import distance

try:
	from muller.clustering.metrics.distance_cache import DistanceCache
	from muller.clustering.methods.changepoints import detect_changepoint
except ModuleNotFoundError:
	from ..metrics.distance_cache import DistanceCache
	from .changepoints import detect_changepoint


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

	@staticmethod
	def _get_cutoff(distances: numpy.ndarray):
		(changepoint_index, changepoint), __ = detect_changepoint(sorted(distances))
		if changepoint_index > len(distances) / 3:
			logger.warning(f"Changepoint not detected.")
			changepoint = pandas.Series([i for i in distances if 0 < i < max(distances)]).quantile(0.05)
		return changepoint


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

	def run(self, pair_array: DistanceCache, starting_genotypes: List[List[str]] = None, similarity_cutoff: Optional[float] = None) -> Tuple[
		List[List[str]], Any]:
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
			if self.cluster_method == 'distance':
				similarity_cutoff = self._get_cutoff(distance_array)
			else:
				similarity_cutoff = 0.05
		logger.debug(f"Using Hierarchical Clustering with similarity cutoff {similarity_cutoff}")

		clusters = self.cluster(reduced_linkage_table, similarity_cutoff, squaremap.index)
		return clusters, linkage_table


def plot_within_cluster_variation():
	from pathlib import Path
	from muller.clustering.metrics import DistanceCalculator
	from muller.clustering.clustercalc import ClusterSet
	import matplotlib.pyplot as plt
	filename = Path("/home/cld100/Documents/github/muller_diagrams/tests/data/tables/real.nature12344-s2.BYB1-G07.xlsx")
	table = pandas.read_excel(filename, sheet_name = 'trajectory').set_index('Trajectory')

	calculator = DistanceCalculator(.03, .97, 'binomial')
	clusterer = HierarchalCluster()

	distances = calculator.run(table)
	x = list()
	y = list()
	for cutoff in range(0, 500):
		logger.info(f"Running cutoff = {cutoff}")
		simcutoff = cutoff / 100
		clusters, _ = clusterer.run(DistanceCache(distances), similarity_cutoff = simcutoff)
		clusters = ClusterSet(clusters, distances)
		k = clusters.calculate_between_cluster_variation()
		logger.debug(f"{simcutoff}, {len(clusters)}, {k}")
		x.append(len(clusters))
		y.append(k)

	plt.scatter(x, y)
	plt.show()


if __name__ == "__main__":
	plot_within_cluster_variation()
