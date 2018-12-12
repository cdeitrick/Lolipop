import argparse
import itertools
from typing import Dict, List, Optional, Tuple

import pandas
from dataclasses import dataclass

try:
	from muller.muller_genotypes.similarity import calculate_pairwise_trajectory_similarity, PairCalculation, PairwiseArrayType
	from muller.muller_genotypes.difference import unlink_unrelated_trajectories
except ModuleNotFoundError:
	from .similarity import calculate_pairwise_trajectory_similarity, PairCalculation, PairwiseArrayType
	from .difference import unlink_unrelated_trajectories

PAIRWISE_P_VALUES: PairwiseArrayType = None
PAIRWISE_CALCULATIONS: Dict[Tuple[str,str], PairCalculation] = None
REMOVED_P_VALUES: PairwiseArrayType = dict()


@dataclass
class GenotypeOptions:
	detection_breakpoint: float  # Minimum frequency to be considered detected.
	fixed_breakpoint: float  # Frequency at which a mutation is considered fixed.
	n_binom: Optional[int]
	similarity_breakpoint: float  # The cutoff indicating two trajectories are related.
	# The cutoff indicating two trajectories that were originally sorted into the same genotype are not
	# actually related.
	difference_breakpoint: float

	@classmethod
	def from_matlab(cls) -> 'GenotypeOptions':
		return GenotypeOptions(
			detection_breakpoint = 0.03,
			fixed_breakpoint = 0.97,
			n_binom = 5,
			similarity_breakpoint = 0.05,
			difference_breakpoint = 0.10
		)

	@classmethod
	def from_breakpoints(cls, detection_cutoff: float, fixed_cutoff: float = None) -> 'GenotypeOptions':
		return GenotypeOptions(
			detection_breakpoint = detection_cutoff,
			fixed_breakpoint = fixed_cutoff if fixed_cutoff else (1 - detection_cutoff),
			n_binom = None,
			similarity_breakpoint = 0.05,
			difference_breakpoint = 0.10
		)

	@classmethod
	def from_parser(cls, parser: argparse.Namespace) -> 'GenotypeOptions':
		compatibility_mode = parser.mode
		detection_breakpoint = float(parser.detection_breakpoint)
		fixed_breakpoint = float(parser.fixed_breakpoint) if parser.fixed_breakpoint else None
		nbinom = int(parser.n_binomial) if hasattr(parser, 'n_binomial') and parser.n_binomial else None
		if fixed_breakpoint is None:
			fixed_breakpoint = 1 - detection_breakpoint
		if compatibility_mode:
			return cls.from_matlab()
		else:
			return GenotypeOptions(
				detection_breakpoint = detection_breakpoint,
				fixed_breakpoint = fixed_breakpoint,
				n_binom = nbinom,
				similarity_breakpoint = parser.similarity_breakpoint,
				difference_breakpoint = parser.difference_breakpoint
			)


@dataclass
class Genotype:
	members: List[str]
	frequency: pandas.Series
	trajectories: pandas.DataFrame

	detected: int  # First timepoint the genotype was detected.
	significant: int  # First significant timepoint
	fixed: Optional[int]  # First fixed timepoint


def _find_genotype_from_trajectory(element: str, all_genotypes: List[List[str]]) -> Optional[List[str]]:
	"""
		Finds the genotype that contains the trajectory.
	Parameters
	----------
	element: int
		The trajectory id.
	all_genotypes: List[List[int]]
		All muller_genotypes that have been calculated.
	Returns
	-------
		The genotype (in the form of a list of trajectory ids) containing the given trajectory id.
		If the trajectory is not contained in any muller_genotypes, returns None
	"""
	candidates = [i for i in all_genotypes if element in i]
	try:
		value = candidates[0]
	except IndexError:
		value = None

	return value


def _group_trajectories_into_genotypes(pairs: PairwiseArrayType, relative_cutoff: float) -> List[List[str]]:
	"""
		Clusters all trajectories into related muller_genotypes.
		By default the first trajectory makes a genotype category
	Parameters
	----------
	pairs: PairwiseArrayType
		A dictionary mapping pairs to p-values.
	relative_cutoff: float; default 0.05
		The cutoff indicating two trajectories are related.

	Returns
	-------
		A list of all muller_genotypes (each of which is a list of ints)
	"""
	if pairs:
		genotype_candidates = [[min(pairs.keys())[0]]]  # by default the first trajectory forms the first genotype.
	else:
		genotype_candidates = []
	seen = set()
	for key, p_value in pairs.items():
		left, right = key
		# ignore pairs that have already been sorted into a genotype.
		if (left, right) in seen or (right, left) in seen:
			continue
		seen.add((left, right))

		if p_value > relative_cutoff:  # are the genotypes related?
			# Check if any of the trajectories are already listed in genotypes.
			# These will return None if no genotype is found.
			genotype_left = _find_genotype_from_trajectory(left, genotype_candidates)
			genotype_right = _find_genotype_from_trajectory(right, genotype_candidates)

			if genotype_left and genotype_right:
				# they are listed under two different genotypes. Combine them.
				if genotype_left != genotype_right:
					genotype_left += genotype_right
					# Remove the redundant genotype.
					genotype_candidates.remove(genotype_right)

			elif genotype_left:
				genotype_left.append(right)
			elif genotype_right:
				genotype_right.append(left)
			else:
				# Neither element is listed. Create a new genotype
				genotype_candidates.append([left, right])
	return genotype_candidates


def calculate_population_genotypes(timeseries: pandas.DataFrame, options: GenotypeOptions) -> List[List[str]]:
	"""
		Clusters trajectories into muller_genotypes.
	Parameters
	----------
	timeseries: pandas.DataFrame
		The timeseries output from import_timeseries.
			- Index -> str
				Names unique to each trajectory.
			- Columns -> int
				The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	options: GenotypeOptions
		The options for this script.
	Returns
	-------
	List[List[str]]
		A list of the muller_genotypes, where each genotype is itself a list of the name of each member trajectory.
		Ex. [
			[A1, A2, A3],
			[B1, B2],
			[C1]
		]
	"""

	# Trajectories represent the population frequencies at each timepoint
	# Each row represents a single timepoint, each column represents a mutation.

	# calculate the similarity between all pairs of trajectories in the population.
	# This is a global variable so that the p-values do not need to be re-computed for every iteration of the genotype filters.
	global PAIRWISE_P_VALUES
	global PAIRWISE_CALCULATIONS
	if PAIRWISE_P_VALUES:
		_current_trajectory_labels = set(timeseries.index)
		pair_array = {k:v for k,v in PAIRWISE_P_VALUES.items() if (k[0] in _current_trajectory_labels and k[1] in _current_trajectory_labels)}

	else:
		pair_array = calculate_pairwise_trajectory_similarity(
			timeseries,
			detection_cutoff = options.detection_breakpoint,
			fixed_cutoff = options.fixed_breakpoint
		)
		if PAIRWISE_CALCULATIONS is None:
			PAIRWISE_CALCULATIONS = {k:v for k,v in pair_array.items() if not isinstance(v, float)}

		pair_array = {k:v.pvalue for k,v in pair_array.items()}

		PAIRWISE_P_VALUES = pair_array

	population_genotypes = _group_trajectories_into_genotypes(pair_array, options.similarity_breakpoint)

	# at the end, look at all trajectories that are not listed and
	# append them as their own category.
	# List of all trajectories
	flattened_genotypes = list(itertools.chain.from_iterable(population_genotypes))
	# Any missing trajectories. Will probably be empty
	other_trajectories = [i for i in timeseries.index if i not in flattened_genotypes]
	population_genotypes.append(other_trajectories)

	# finally, for each genotype, make sure each trajectory pair has some
	# non-trivial linkage (say, >0.0005). this avoids falsely linking together
	# trajectories. if not, divide the offending trajectories into two
	# camps, and sort the rest according to which one they are more
	# closely linked to. repeat until everything is linked.
	# pprint(pair_array)

	while True:
		starting_size_of_the_genotype_array = len(population_genotypes)
		population_genotypes = unlink_unrelated_trajectories(population_genotypes[:], pair_array, options.difference_breakpoint)
		if len(population_genotypes) == starting_size_of_the_genotype_array:
			break

	# all_genotypes[population_id] = population_genotypes
	# TODO Fix so that it returns a genotype for each population
	return [i for i in population_genotypes if i]  # Only return non-empty lists.


def _calculate_mean_frequencies_of_trajectories(name: str, genotype_timeseries: pandas.DataFrame, genotype: List[str]) -> pandas.Series:
	mean_genotype_timeseries = genotype_timeseries.mean()
	mean_genotype_timeseries['members'] = "|".join(map(str, genotype))
	mean_genotype_timeseries.name = name

	return mean_genotype_timeseries


def calculate_mean_genotype(all_genotypes: List[List[str]], timeseries: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Calculates the mean frequency of each genotype ate every timepoint.
	Parameters
	----------
	all_genotypes: List[Genotype]
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
		genotype_timeseries = timeseries.loc[genotype]
		mean_genotype_timeseries = _calculate_mean_frequencies_of_trajectories(f"genotype-{index}", genotype_timeseries, genotype)
		mean_genotypes.append(mean_genotype_timeseries)

	mean_genotypes = pandas.DataFrame(mean_genotypes)
	return mean_genotypes


def workflow(timepoints: pandas.DataFrame, options: GenotypeOptions) -> pandas.DataFrame:
	"""

	Parameters
	----------
	timepoints: pandas.DataFrame
		A timeseries dataframe, usually generated from `import_table.import_trajectory_table`.
			- Index -> str
				Names unique to each trajectory.
			- Columns -> int
				The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	options: GenotypeOptions
		An instance of `GenotypeOptions` with the desired options. Can be generated automatically via `GenotypeOptions.from_breakpoints(0.03)`.

	Returns
	-------
	pandas.DataFrame
	"""

	genotypes = calculate_population_genotypes(timepoints, options)

	_mean_genotypes = calculate_mean_genotype(genotypes, timepoints)

	# _mean_genotypes.to_csv(str(filename.with_suffix('.mean.tsv')), sep = '\t')
	return _mean_genotypes


if __name__ == "__main__":
	pass
