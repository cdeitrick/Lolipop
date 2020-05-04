""" Implements a class which can compute statistical values over a set of clusters. Implemented as its own class to make is usable elsewhere."""
import itertools
import statistics
from typing import Dict, List, Tuple, Union

from loguru import logger

from muller.clustering.metrics import DistanceCache

ClusterType = List[str]


class ClusterSet:
	""" Provides descriptive statistics for a set of clusters."""
	def __init__(self, clusters: List[ClusterType], distances: Union[DistanceCache, Dict[Tuple[str, str], float]]):
		if not isinstance(distances, DistanceCache):
			self.distances = DistanceCache(distances)
		else:
			self.distances = distances

		self.number_of_points = len(self.distances.squareform().index)

		self.clusters = clusters
		self._genotypes = None  # Cache for the `self.genotypes` property.

	def __len__(self) -> int:
		return len(self.clusters)

	def __getitem__(self, item):
		return self.genotypes[item]

	def calculate_index(self) -> float:
		""" Calculates the CH index for the clusters."""

		B = self.calculate_within_cluster_variations()
		W = self.calculate_between_cluster_variation()
		logger.debug(f"{B}, {W}")
		result = (B / len(self) - 1) / (W / (self.number_of_points - len(self)))
		return result

	def calculate_between_cluster_variation(self):
		clusters = self.get_genotypes(self.clusters)
		combinations = itertools.combinations(clusters.keys(), 2)
		cluster_distances = list()
		for left_genotype, right_genotype in combinations:
			# Need to calculate the distance between these two clusters
			left_cluster = clusters[left_genotype]
			right_cluster = clusters[right_genotype]

			distance = self.calculate_distance_between_clusters(left_cluster, right_cluster, 'min')
			#logger.debug(f"{left_genotype}, {right_genotype}, ({len(left_cluster), len(right_cluster)}), {distance}")
			cluster_distances.append(distance)

		mean_cluster_distance = statistics.mean(cluster_distances)
		variances = [(i - mean_cluster_distance) ** 2 for i in cluster_distances]
		return sum(variances) / len(clusters)

	def calculate_distance_between_clusters(self, left_cluster: ClusterType, right_cluster: ClusterType, method = 'min') -> float:
		# Need to pair the left trajectories with the right trajectories, but not left to left.

		between_cluster_distances = [self.distances[left, right] for right in right_cluster for left in left_cluster]

		if method == 'min':
			return min(between_cluster_distances)
		elif method == 'max':
			return max(between_cluster_distances)
		elif method == 'mean':
			return statistics.mean(between_cluster_distances)
		else:
			message = f"Invalid method for calculating the between-cluster distance: '{method}'"
			raise ValueError(message)

	def calculate_silhouette_coefficient(self, label: str, cluster_label: str, other_label: str):
		""" Calculates the Silhouette Coefficient for a single point.
			Parameters
			----------
			label: str
				THe label of the sample
			cluster_label: str
				The label of the cluster that `label` belongs to.
			other_label: str
				The label of the next closest cluster.
		"""
		pass

	def calculate_silhouette_coefficients(self) -> Dict[str, float]:
		coefficients = dict()
		for cluster_name, cluster in self.genotypes.items():
			if len(cluster) == 1: continue
			closest_cluster = self.get_closest_cluster(cluster_name)
			closest_cluster_members = self[closest_cluster]
			for member in cluster:
				distances_between_other_members = [self.distances[member, other] for other in cluster]
				distances_between_closest_cluster = [self.distances[member, other] for other in closest_cluster_members]
				# -1 sincec `member` is included as 0
				average_distance = sum(distances_between_other_members) / (len(distances_between_other_members) - 1)
				average_distance_closest = sum(distances_between_closest_cluster) / len(distances_between_closest_cluster)
				c = (average_distance_closest - average_distance) / (max(distances_between_closest_cluster + distances_between_other_members))
				coefficients[member] = c
		return coefficients

	@property
	def genotypes(self) -> Dict[str, ClusterType]:
		if self._genotypes is None:
			self._genotypes = self.get_genotypes(self.clusters)
		return self._genotypes

	@staticmethod
	def calculate_within_cluster_variation(cluster: ClusterType, pair_array: DistanceCache) -> float:
		""" Row labels and column names are samples."""

		if len(cluster) > 1:
			cluster_distances = [pair_array[i, j] for i, j in itertools.combinations(cluster, 2)]
		else:
			cluster_distances = [0]
		cluster_mean = statistics.mean(cluster_distances)

		variances = [(i - cluster_mean) ** 2 for i in cluster_distances]
		W = sum(variances) / len(cluster)
		return W

	def calculate_within_cluster_variations(self) -> float:
		""" Calulates the total within-cluster variation for all clusters in `clusters`.
		"""
		cluster_variations = [self.calculate_within_cluster_variation(cluster, self.distances) for cluster in self.clusters]
		return sum(cluster_variations)

	def get_genotypes(self, clusters: List[ClusterType] = None) -> Dict[str, ClusterType]:
		""" Generates labels for each cluster."""
		if clusters is None: clusters = self.clusters
		clusterdict = dict()
		for index, cluster in enumerate(clusters, start = 1):
			# Test if each trajectory shares a common name. Usually due to testing genotypes.
			try:
				_names = {i.split('-')[1] for i in cluster}
			except IndexError:
				_names = set()
			if len(_names) == 1:
				genotype_name = "genotype-" + list(_names)[0]
			else:
				genotype_name = f"genotype-{index}"
			clusterdict[genotype_name] = cluster
		return clusterdict

	def get_closest_cluster(self, label: str) -> str:
		""" Return the label of the closest cluster to `label`"""
		dist = dict()
		for name, cluster in self.genotypes.items():
			if name == label: continue
			d = self.calculate_distance_between_clusters(self[label], self[name])
			dist[name] = d
		return min(dist.items(), key = lambda s: s[1])[0] 
