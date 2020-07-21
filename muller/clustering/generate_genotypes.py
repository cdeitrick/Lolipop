import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger

try:
	from muller.clustering import metrics, genotype_reorder, hierarchy
	from muller import filters
	from muller.dataio import projectdata
except ModuleNotFoundError as exception:
	logger.warning(str(exception))
	from .. import filters
	from . import metrics, hierarchy
def is_trajectory_labeled_by_genotype(label):
	regex = "trajectory-[a-z]+-[0-9]+"
	match = re.search(regex, label)
	return match

def generate_genotype_name(index: int, members: List[str]):
	""" Generates the genotype names.
		When the trajectories are named like "trajectory-aqua-2", the genotype will be named "genotype-aqua".
	"""
	# Check if the trajectories follow the naming convention
	trajectory_label_results = [is_trajectory_labeled_by_genotype(m) for m in members]
	all_trajectories_follow_convention = all(trajectory_label_results)
	try:
		unique_genotype_labels = sorted(set(i.split('-')[1] for i in members))
	except IndexError:
		unique_genotype_labels = list()
	all_trajectories_specify_same_genotype = len(unique_genotype_labels) == 1

	if all_trajectories_follow_convention and all_trajectories_specify_same_genotype:
		genotype_name = unique_genotype_labels[0]
	else:
		genotype_name = str(index)
	genotype_label = f"genotype-{genotype_name}"
	return genotype_label

class ClusterMutations:
	"""
	Parameters
	----------
	dlimit, flimit:float
		The detection, and fixed cutoff values, respectively
	metric: {'binomial', 'pearson', 'minkowski'}
		The distance metric to determine trajectory similarity.
	starting_genotypes: List[List[str]]
		A list of genotypes to start the clustering algorithm with. The distance metrics will be modified so that the trajectories specified are
		grouped together.
	"""

	def __init__(self, metric: str, dlimit: float, flimit: float,
			starting_genotypes: Optional[List[List[str]]] = None, threads: Optional[int] = None):
		self.metric: str = metric
		self.dlimit: float = dlimit
		self.flimit: float = flimit
		self.known_genotypes: List[List[str]] = starting_genotypes if starting_genotypes else []
		self.filename_pairwise = None
		self.pairwise_distances_full = None # overwritten in self.get_pairwise_distances.

		#self.filename_pairwise = filename_pairwise # Could be used to reuse a table of pairwise distances.
		# The `breakpoints` value is a bit arbitrary, so it should be safe to hard-code it.
		# 	This will actually prevent the most common error when sorting genotypes (i.e. no breakpoints given) so it's worth
		#	hard-coding it to prevent that issue.
		self.breakpoints = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

		self.distance_calculator = metrics.DistanceCalculator(
			detection_limit = self.dlimit,
			fixed_limit = self.flimit,
			metric = self.metric,
			threads = threads
		)

		self.clusterer = hierarchy.HierarchalCluster()

		self.organizer = genotype_reorder.SortGenotypeTableWorkflow(
			dlimit = dlimit,
			flimit = flimit
		)

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
		mean_genotype_timeseries.index = [int(i) for i in mean_genotype_timeseries.index]
		mean_genotype_timeseries['members'] = "|".join(map(str, genotype))
		mean_genotype_timeseries.name = name

		return mean_genotype_timeseries

	@staticmethod
	def _load_pairwise_distances(filename: Path) -> Dict[Tuple[str, str], float]:
		""" Reads pre-computed pairwise distances from a previous run. Typically found in the /tables/.distance.tsv table."""
		table_distance_pairwise = pandas.read_csv(filename, sep = "\t")

		pair_array = dict()
		for index, row in table_distance_pairwise.iterrows():
			values = {(index, i): row[i] for i in table_distance_pairwise.columns}
			values_inverse = {k[::-1]: v for k, v in
				values.items()}  # Faster to include reversed keys rather than trying to get (l,r) and (r,l) keys.

			pair_array.update(values)
			pair_array.update(values_inverse)

		return pair_array




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
			#genotype_name = generate_genotype_name(index, genotype_timeseries.index)
			genotype_name = f"genotype-{index}"
			mean_genotype_timeseries = self._calculate_mean_frequencies_of_trajectories(genotype_name, genotype_timeseries, genotype)
			mean_genotypes.append(mean_genotype_timeseries)

		mean_genotypes = pandas.DataFrame(mean_genotypes)
		# For consistency
		mean_genotypes.index.name = 'Genotype'

		# Place the `members` column in the leftmost column
		mean_genotypes = mean_genotypes[['members'] + [i for i in mean_genotypes if i != 'members']]

		return mean_genotypes

	def generate_genotype_table(self, timepoints: pandas.DataFrame, genotypes: List[List[str]]) -> Tuple[pandas.DataFrame, Dict[str,List[str]]]:

		mean_genotypes = self.calculate_mean_genotype(genotypes, timepoints)

		# Try to keep the genotype members separate so that `mean_genotypes` is consistent.
		# The table should only have the timeseries frequency values and be indexed by genotype label.
		genotype_members = mean_genotypes.pop('members')
		# Convert the `genotype_members` pandas.Series object to a dictionary mapping genotypes to member trajectories.
		genotype_members = {k: v.split('|') for k, v in genotype_members.items()}

		return mean_genotypes, genotype_members

	def get_pairwise_distances(self, trajectories: pandas.DataFrame):
		if self.filename_pairwise:
			pair_array = self._load_pairwise_distances(self.filename_pairwise)
		else:
			pair_array = self.distance_calculator.run(trajectories)

		self.pairwise_distances_full = metrics.DistanceCache(pair_array)  # Keep a record of the pairwise distances before filtering.
		return self.pairwise_distances_full

	def run(self, trajectories: pandas.DataFrame, distance_cutoff:Optional[float] = None) -> projectdata.DataGenotypeInference:
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
		distance_cutoff: Optional[float]
			Used to determine the distance cutoff when clustering genomtypes.
		"""

		modified_trajectories = trajectories.copy(deep = True)  # To avoid unintended changes

		# Calculate the pairwise distances between each pair of mutational trajectories.
		pairwise_distances = self.get_pairwise_distances(modified_trajectories)

		# Calculate the genotypes
		cluster_result = self.clusterer.run(
			pairwise_distances,
			starting_genotypes = self.known_genotypes,
			similarity_cutoff = distance_cutoff
		)
		genotype_table, genotype_members = self.generate_genotype_table(modified_trajectories, cluster_result.clusters)

		sorted_genotype_table = self.organizer.run(genotype_table)
		output_data = projectdata.DataGenotypeInference(
			table_trajectories = modified_trajectories,
			table_genotypes = sorted_genotype_table,
			genotype_members = genotype_members,
			matrix_distance = pairwise_distances,
			clusterdata = cluster_result,
			table_trajectories_info = None
		)
		print(genotype_table.to_string())
		return output_data

