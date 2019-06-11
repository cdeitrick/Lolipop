from typing import List, Optional, Tuple

import numpy
import pandas
from loguru import logger
try:
	from muller.clustering import filters, metrics, methods
except ModuleNotFoundError:
	from . import filters, metrics, methods


class ClusterMutations:
	"""
	Parameters
	----------
	dlimit, flimit:float
		The detection, significant, and fixed cutoff values, respectively
	sbreakpoint: float
		The cutoff value to use when clustering trajectories into genotypes. Two trajectories must have a distance less than this value to be
		considered members of the same genotype.
	dbreakpoint: float
		Only used when using the two-step method. Governs when the genotype contains two trajectories which should be split into separate genotypes.
	method: {'two-step', 'hierarchy'}
		The clustering method to use. The hierarchical method uses the maximum distance between a point and a candidate cluster to determine whether
		to group said point to the cluster.
	metric: {'binomial', 'pearson', 'minkowski'}
		The distance metric to determine trajectory similarity.
	starting_genotypes: List[List[str]]
		A list of genotypes to start the clustering algorithm with. The distance metrics will be modified so that the trajectories specified are
		grouped together.
	Returns
	-------
	pandas.DataFrame, pandas.Series, numpy.array
		- The genotype table
		- A map of genotypes to members.
		- A linkage matrix, if hierarchical clustering was used.
	"""

	def __init__(self, method: str, metric: str, dlimit: float, flimit: float, sbreakpoint: float, dbreakpoint: float, breakpoints: List[float],
			starting_genotypes: List[List[str]], include_single:bool):
		self.method: str = method
		self.metric: str = metric
		self.dlimit: float = dlimit
		self.flimit: float = flimit
		self.sbreakpoint: float = sbreakpoint
		self.dbreakpoint: float = dbreakpoint
		self.breakpoints: List[float] = breakpoints
		self.starting_genotypes: List[List[str]] = starting_genotypes
		self.include_single = include_single

		# These will be updated when run() is called.
		self.pairwise_distances: metrics.PairwiseCalculationCache = metrics.PairwiseCalculationCache()  # Empty cache that will be replaced in generate_genotypes().
		self.genotype_table = self.genotype_members = self.linkage_table = self.rejected_trajectories = None

	def run(self, trajectories: pandas.DataFrame):
		"""
			Run the genotype clustering workflow.
		Parameters
		----------
		trajectories: pandas.DataFrame
			A timeseries dataframe, usually generated from `import_table.import_trajectory_table`.
				- Index -> str
					Names unique to each trajectory.
				- Columns -> int
					The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
					Each trajectory/timepoint will include the observed frequency at each timepoint.
		"""
		trajectory_filter = filters.TrajectoryFilter(detection_cutoff = self.dlimit, fixed_cutoff = self.flimit, exclude_single = not self.include_single)
		genotype_filter = filters.GenotypeFilter(detection_cutoff = self.dlimit, fixed_cutoff = self.flimit, frequencies = self.breakpoints)
		modified_trajectories = trajectories.copy(deep = True)  # To avoid unintended changes
		if self.breakpoints:
			# The filters should run
			_iterations = 20  # arbitrary, used to make sure the program does not encounter an infinite loop.
			modified_trajectories = trajectory_filter.run(modified_trajectories)
		else:
			# The filters are disabled.
			_iterations = 0  # The for loop shouldn't excecute.

		# Calculate the distance between all possible pair of trajectories.
		pair_array = metrics.calculate_pairwise_metric(
			modified_trajectories,
			detection_cutoff = self.dlimit,
			fixed_cutoff = self.flimit,
			metric = self.metric
		)
		self.pairwise_distances = metrics.PairwiseCalculationCache(pair_array)

		# Calculate the initial genotypes
		genotype_table, genotype_members, linkage_matrix = self.generate_genotypes(modified_trajectories)

		rejected_members = {i:'filtered-genotype-0' for i in trajectories.index if i not in modified_trajectories}
		for index in range(_iterations):
			invalid_members = genotype_filter.run(genotype_table, genotype_members)
			if invalid_members:
				for i in invalid_members:
					rejected_members[i] = f"filtered-genotype-{index+1}"
				# Remove these trajectories from the trajectories table.
				modified_trajectories = modified_trajectories[~modified_trajectories.index.isin(invalid_members)]

				# Need to remove trajectories from the distance matrix so they are not included in the clustering method.
				self.pairwise_distances.reduce(modified_trajectories.index)
				# Re-calculate the genotypes based on the remaining trajectories.

				genotype_table, genotype_members, linkage_matrix = self.generate_genotypes(modified_trajectories)
			else:
				break
		# Build a table of trajectories that were rejected.
		self.rejected_trajectories = trajectories.loc[sorted(rejected_members.keys())]
		self.rejected_trajectories['genotype'] = [rejected_members[i] for i in self.rejected_trajectories.index]

		self.genotype_table = genotype_table
		self.genotype_members = genotype_members
		self.linkage_table = linkage_matrix

		return genotype_table, genotype_members

	def generate_genotypes(self, timepoints: pandas.DataFrame) -> Tuple[pandas.DataFrame, pandas.Series, Optional[numpy.array]]:
		if self.method == "matlab" or self.method == 'twostep':
			genotypes = methods.twostep_method(timepoints, self.pairwise_distances, self.sbreakpoint, self.dbreakpoint, self.starting_genotypes)
			linkage_matrix = None
		elif self.method == "hierarchy":
			genotypes, linkage_matrix = methods.hierarchical_method(self.pairwise_distances, self.sbreakpoint, starting_genotypes = self.starting_genotypes)
		else:
			raise ValueError(f"Invalid clustering method: {self.method}")

		mean_genotypes = self.calculate_mean_genotype(genotypes, timepoints)
		genotype_members = mean_genotypes.pop('members')
		return mean_genotypes, genotype_members, linkage_matrix

	@staticmethod
	def _calculate_mean_frequencies_of_trajectories(name: str, genotype_timeseries: pandas.DataFrame, genotype: List[str]) -> pandas.Series:
		"""
			Generates a mean timeseries for a genotype given the timeseries of member trajectories.
		Parameters
		----------
		name: str
			The name to assign to this genotype. Usually just 'genotype-' forllowed by an integer.
		genotype_timeseries: pandas.DataFrame
			A dataframe of member trajectory timeseries.
		genotype: List[str]
			A list of the labels of all member trajectories. Used to track the which trajectories belong to this genotype.

		Returns
		-------
		pandas.DataFrame
			The calculated mean of all member trajectories at each timepoint.
		"""
		mean_genotype_timeseries = genotype_timeseries.mean()
		mean_genotype_timeseries['members'] = "|".join(map(str, genotype))
		mean_genotype_timeseries.name = name

		return mean_genotype_timeseries

	def calculate_mean_genotype(self, all_genotypes: List[List[str]], timeseries: pandas.DataFrame) -> pandas.DataFrame:
		"""
			Calculates the mean frequency of each genotype ate every timepoint.
		Parameters
		----------
		all_genotypes: List[List[str]
			A list of all muller_genotypes for a given population
		timeseries: pandas.DataFrame
			Must have a 'Trajectory' column along with the columns of the original table that represent timepoints.
			Each row corresponds to a single mutational trajectory.

		Returns
		-------
		A dataframe where every row corresponds to a genotype.
		member trajectories are listed under the 'members' column.
		every column represents a timepoint.
		"""
		mean_genotypes = list()
		for index, genotype in enumerate(all_genotypes, start = 1):
			try:
				genotype_timeseries = timeseries.loc[genotype]
			except Exception as exception:
				logger.critical(f"Missing Trajectory Labels: {genotype} - {timeseries.index}")
				raise exception
			mean_genotype_timeseries = self._calculate_mean_frequencies_of_trajectories(f"genotype-{index}", genotype_timeseries, genotype)
			mean_genotypes.append(mean_genotype_timeseries)

		mean_genotypes = pandas.DataFrame(mean_genotypes)
		# For consistency
		mean_genotypes.index.name = 'Genotype'

		return mean_genotypes
