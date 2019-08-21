import itertools
from typing import Dict, List, Optional, Tuple
from pathlib import Path
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
	filename_pairwise: Optional[Path]
		Path to a tables/.distance.tsv file from a previous run using the same input parameters and table.
	Returns
	-------
	pandas.DataFrame, pandas.Series, numpy.array
		- The genotype table
		- A map of genotypes to members.
		- A linkage matrix, if hierarchical clustering was used.
	"""

	def __init__(self, method: str, metric: str, dlimit: float, flimit: float, sbreakpoint: float, dbreakpoint: float, breakpoints: List[float],
			starting_genotypes: List[List[str]], trajectory_filter: filters.TrajectoryFilter, filename_pairwise:Optional[Path] = None, threads:Optional[int] = None):
		self.method: str = method
		self.metric: str = metric
		self.dlimit: float = dlimit
		self.flimit: float = flimit
		self.sbreakpoint: float = sbreakpoint
		self.dbreakpoint: float = dbreakpoint
		self.breakpoints: List[float] = breakpoints
		self.starting_genotypes: List[List[str]] = starting_genotypes
		self.trajectory_filter = trajectory_filter
		self.filename_pairwise = filename_pairwise
		# These will be updated when run() is called.
		self.pairwise_distances: metrics.DistanceCache = metrics.DistanceCache()  # Empty cache that will be replaced in generate_genotypes().
		self.pairwise_distances_full: metrics.DistanceCache = metrics.DistanceCache()  # Empty cache that will be replaced in generate_genotypes().
		self.genotype_table = self.genotype_members = self.linkage_table = self.rejected_trajectories = None

		self.distance_calculator = metrics.DistanceCalculator(self.dlimit, self.flimit, self.metric, threads = threads)
		self.genotype_filter = filters.GenotypeFilter(
			detection_cutoff = self.dlimit,
			fixed_cutoff = self.flimit,
			frequencies = self.breakpoints
		)
	@staticmethod
	def _load_pairwise_distances(filename:Path)->Dict[Tuple[str,str], float]:
		""" Reads pre-computed pairwise distances from a previous run. Typically found in the /tables/.distance.tsv table."""
		table_distance_pairwise = pandas.read_csv(filename, sep = "\t")

		pair_array = dict()
		for index, row in table_distance_pairwise.iterrows():
			values = {(index,i): row[i] for i in table_distance_pairwise.columns}
			values_inverse = {k[::-1]:v for k,v in values.items()} # Faster to include reversed keys rather than trying to get (l,r) and (r,l) keys.

			pair_array.update(values)
			pair_array.update(values_inverse)

		return pair_array

	def get_pairwise_distances(self, trajectories:pandas.DataFrame):

		if self.filename_pairwise:
			pair_array = self._load_pairwise_distances(self.filename_pairwise)
		else:
			pair_array = self.distance_calculator.run(trajectories)


		self.pairwise_distances_full = metrics.DistanceCache(pair_array) # Keep a record of the pairwise distances before filtering.
		return self.pairwise_distances_full

	def run(self, trajectories: pandas.DataFrame) -> Tuple[pandas.DataFrame, Dict[str, List[str]]]:
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

		modified_trajectories = trajectories.copy(deep = True)  # To avoid unintended changes
		if self.trajectory_filter:
			modified_trajectories = self.trajectory_filter.run(modified_trajectories)
		if self.breakpoints:
			# The filters should run
			_iterations = 20  # arbitrary, used to make sure the program does not encounter an infinite loop.
		else:
			# The filters are disabled.
			_iterations = 0  # The for loop shouldn't excecute.

		# Calculate the distance between all possible pair of trajectories.

		self.pairwise_distances =self.get_pairwise_distances(modified_trajectories)

		# Calculate the initial genotypes
		genotype_table, genotype_members, linkage_matrix = self.generate_genotypes(modified_trajectories)

		for index in range(_iterations):
			invalid_members = self.genotype_filter.run(genotype_table, genotype_members)
			if invalid_members:
				# Remove these trajectories from the trajectories table.
				modified_trajectories = modified_trajectories[~modified_trajectories.index.isin(invalid_members)]

				# Need to remove trajectories from the distance matrix so they are not included in the clustering method.
				self.pairwise_distances.reduce(modified_trajectories.index)
				# Re-calculate the genotypes based on the remaining trajectories.

				genotype_table, genotype_members, linkage_matrix = self.generate_genotypes(modified_trajectories)
			else:
				break

		self.genotype_table = genotype_table
		self.linkage_table = linkage_matrix

		self.genotype_members = {k: v.split('|') for k, v in genotype_members.items()}

		nonfiltered_trajectories = list(itertools.chain.from_iterable(self.genotype_members.values()))
		filtered_trajectories = [i for i in trajectories.index if i not in nonfiltered_trajectories]

		self.genotype_members['genotype-filtered'] = filtered_trajectories

		return self.genotype_table, self.genotype_members

	def generate_genotypes(self, timepoints: pandas.DataFrame) -> Tuple[pandas.DataFrame, pandas.Series, Optional[numpy.array]]:
		if self.method == "matlab" or self.method == 'twostep':
			genotypes = methods.twostep_method(timepoints, self.pairwise_distances, self.sbreakpoint, self.dbreakpoint, self.starting_genotypes)
			linkage_matrix = None
		elif self.method == "hierarchy":
			genotypes, linkage_matrix = methods.hierarchical_method(self.pairwise_distances, self.sbreakpoint,
				starting_genotypes = self.starting_genotypes)
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
