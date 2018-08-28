import pandas

import itertools
import math
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple, Union
from pathlib import Path
import argparse

try:
	from muller.time_series_import import import_timeseries
except ModuleNotFoundError:
	from time_series_import import import_timeseries

# The data structure of a agenotype
Genotype = List[int]


@dataclass
class PairArrayValue:
	"""
		Holds the values assigned to parray in the original script.
	"""
	population: int
	left: int
	right: int
	p_value: float
	sigma_value: float
	difbar: float


PairwiseArrayType = Dict[Tuple[int, int], PairArrayValue]


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
	def from_parser(cls, parser:argparse.Namespace) -> 'GenotypeOptions':
		compatibility_mode = parser.mode
		detection_breakpoint = float(parser.detection_breakpoint)
		fixed_breakpoint = float(parser.fixed_breakpoint) if parser.fixed_breakpoint else None
		nbinom = int(parser.n_binomial) if hasattr(parser, 'n_binomial') and parser.n_binomial else None
		if fixed_breakpoint is None:
			fixed_breakpoint = 1-detection_breakpoint
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

def calculate_pairwise_similarity(population_id: int, timeseries: pandas.DataFrame, detection_cutoff: float,
		fixed_cutoff: float, n_binomial: int) -> PairwiseArrayType:
	"""
		Calculates the relative similarity between all trajectory pairs.

		these are consistent with n independent draws from a normal distribution,
        assuming an average variance of  n_{binom} p (1-p), where p is the
        average of the two frequencies, and n_{binom} is picked arbitrarily as a
        value that gives reasonable uncertainties given our data

        n random draws will have
        \sigma_{tot}^2 = \sum_i \sigma_i^2 =  n_{binom} \sum_i p_i(1 - p_i), from
        a property of normal distributions

        for \bar{X} = \sum_{i}^{n_X} X/n_X

        \sigma_{\bar{X}} = \sigma_{tot}/n_X = \sqrt( \sum_i p_i(1-p_i))/\sqrt{n_X}

        Finally, given \sigma and \bar{X} we can construct a p-value for
        the measurement by numerically computing the following integral:

        1 - \int_{- \bar{X}/ \sigma_{\bar{X}}}^{\bar{X}/
        \sigma_{\bar{X}}} \frac{1}{\sqrt{2 \pi}} e^{-x^2/2} dx


	Parameters
	----------
	population_id
	timeseries: pandas.DataFrame
	detection_cutoff: float
	fixed_cutoff: float
	n_binomial: int

	Returns
	-------
	dict of PairArrayValue
	Each key in the dictionary corresponds to a pair of trajectory ids which map to the p-value for that pair.
	The order of ids does not matter.

	"""

	# get a list of all possible trajectory pairs, ignoring element order.
	trajectories = timeseries[[i for i in timeseries.columns if i not in ['Population', 'Position', 'Trajectory']]]
	combos = sorted(itertools.combinations(timeseries['Trajectory'].values, 2))
	pair_array = dict()

	for pair in combos:
		# Unpack the trajectory labels/ids
		left, right = pair

		# Extract the trajectory frequencies for both trajectories in the pair.
		# trajectory indicies start at 1, need to convert to 0-indexed.
		left_trajectories = trajectories.iloc[left - 1]
		right_trajectories = trajectories.iloc[right - 1]

		# Merge into a dataframe for convienience
		df = pandas.concat([left_trajectories, right_trajectories], axis = 1)

		# Calculate similarity
		# Remove timepoints where at least one trajectory was not fixed or undetected.
		filtered_df = df[(df < fixed_cutoff).any(axis = 1)]
		filtered_df = filtered_df[(filtered_df > detection_cutoff).any(axis = 1)]

		if filtered_df.empty:
			# Both trajectories have no timepoints where they are detected but not yet fixed.
			# Assign a p-value of 0.0
			pair_array_value = PairArrayValue(
				population_id, pair[0], pair[1], math.nan, math.nan, math.nan
			)
		else:
			# Find the mean frequency of each timepoint
			# index is timepoints,  values are frequencies
			mean_frequencies: pandas.Series = filtered_df.mean(axis = 1)

			# Calculate sigma_total

			n_binom = n_binomial if n_binomial else len(mean_frequencies)

			# noinspection PyTypeChecker
			sigma_freq: pandas.Series = (mean_frequencies * (
					1 - mean_frequencies)) / n_binom  # The 5 is arbitrary see above

			# Difference of frequencies at each timepoint
			difference: pandas.Series = filtered_df.iloc[:, 0] - filtered_df.iloc[:, 1]

			sigmapair: float = math.sqrt(sigma_freq.sum()) / len(difference)

			# Sum of differences
			difbar: float = abs(difference).sum() / len(difference)

			# Calculates the relative similarity between all trajectory pairs.
			X = difbar / (math.sqrt(2) * sigmapair)
			pval: float = 1 - math.erf(X)

			pair_array_value = PairArrayValue(
				population_id, pair[0], pair[1], pval, sigmapair, difbar
			)
		# Add the p-value information for the forward and reverse pairs (will be the same, no need to do again)
		pair_array[(left, right)] = pair_array_value
		pair_array[(right, left)] = pair_array_value
	return pair_array


def find_genotype(element: int, all_genotypes: List[List[int]]) -> Optional[Genotype]:
	"""
		Finds the genotype that contains the trajectory.
	Parameters
	----------
	element: int
		The trajectory id.
	all_genotypes: List[List[int]]
		All genotypes that have been calculated.
	Returns
	-------
		The genotype (in the form of a list of trajectory ids) containing the given trajectory id.
		If the trajectory is not contained in any genotypes, returns None
	"""
	candidates = [i for i in all_genotypes if element in i]
	try:
		value = candidates[0]
	except IndexError:
		value = None

	return value


def find_all_genotypes(pairs: PairwiseArrayType, relative_cutoff: float) -> List[Genotype]:
	"""
		Clusters all trajectories into related genotypes.
		By default the first trajectory makes a genotype category
	Parameters
	----------
	pairs: PairArrayValue
		A dictionary mapping pairs to p-values.
	relative_cutoff: float; default 0.05
		The cutoff indicating two trajectories are related.

	Returns
	-------
		A list of all genotypes (each of which is a list of ints)
	"""
	genotype_candidates = [[1]]  # by default the first trajectory forms the first genotype.
	seen = set()
	for pair in pairs.values():
		# ignore pairs that have already been sorted into a genotype.
		if (pair.left, pair.right) in seen or (pair.right, pair.left) in seen:
			continue
		seen.add((pair.left, pair.right))

		if pair.p_value > relative_cutoff:  # are the genotypes related?
			# Check if any of the trajectories are already listed in genotypes.
			# These will return None if no genotype is found.
			genotype_left = find_genotype(pair.left, genotype_candidates)
			genotype_right = find_genotype(pair.right, genotype_candidates)

			if genotype_left and genotype_right:
				# they are listed under two different genotypes. Combine them.
				if genotype_left != genotype_right:
					genotype_left += genotype_right
					# Remove the redundant genotypes.
					genotype_candidates.remove(genotype_right)

			elif genotype_left:
				genotype_left.append(pair.right)
			elif genotype_right:
				genotype_right.append(pair.left)
			else:
				# Neither element is listed. Create a new genotype
				genotype_candidates.append([pair.left, pair.right])
	return genotype_candidates


def get_p_value(left: int, right: int, pairwise_array: PairwiseArrayType) -> Optional[PairArrayValue]:
	"""	Finds the p-value for a given pair."""
	return pairwise_array[(left, right)]


def split_genotype_in_two(genotype: Genotype, unlinked_trajectories: pandas.DataFrame,
		pair_array: PairwiseArrayType, link_cut: float) -> Tuple[Genotype, Genotype]:
	"""
		Splits a genotype into smaller genotypes if some members are not related to some other members. This may happen
		when a member is included into a genotype due to its paired member but ends up not related to some of
		the other members that where found later.
	Parameters
	----------
	genotype: List[int]
		The genotype to split.
	unlinked_trajectories: pandas.DataFrame
		A dataframe of trajectory pairs and their corresponding p-value.
		Columns:
			- left
			- right
			- pvalue
	pair_array
		The list of p-values for each pair.
	link_cut: float
		The cuttoff value to choose whether a member is related to the other members.

	Returns
	-------
		a tuple of two genotypes.
	"""
	# Find the index of the minimum p-value after subtracting the link cutoff.
	minimum_pvalue_index = (unlinked_trajectories['pvalue'] - link_cut).abs().idxmin()
	# Get the row with the identified minimu mp-value.
	minimum_pvalue_genotype = unlinked_trajectories.loc[minimum_pvalue_index]

	# Form two new genotypes based on the two trajectories corresponding to the minimum p-value
	new_genotype_1_base = int(minimum_pvalue_genotype['left'])
	new_genotype_2_base = int(minimum_pvalue_genotype['right'])

	# Make sure genotype 1 includes the lower id-value. Not important, but maintains parity with matlab script.
	new_genotype_1_base, new_genotype_2_base = sorted([new_genotype_1_base, new_genotype_2_base])
	new_genotype_1 = [new_genotype_1_base]
	new_genotype_2 = [new_genotype_2_base]

	# Sort each member in the current genotype into one of the new genotypes.
	for genotype_member in genotype:
		genotype_member = int(genotype_member)
		# Check if the current genotype member is already contained in one of the genotypes.
		# Should only be one of the two trajectories used to form a new genotype.
		if genotype_member in new_genotype_1 or genotype_member in new_genotype_2:
			pass
		else:
			# Use the highest p-value to determine which genotype to add the member to.
			# P-values should correspond to the current member and the base member of the new genotypes.
			p_value_1 = get_p_value(new_genotype_1_base, genotype_member, pair_array)
			p_value_2 = get_p_value(new_genotype_2_base, genotype_member, pair_array)

			if p_value_1.p_value >= p_value_2.p_value:
				new_genotype_1.append(genotype_member)
			else:
				new_genotype_2.append(genotype_member)
	return new_genotype_1, new_genotype_2


def split_unlinked_genotypes(all_genotypes: List[List[int]], pair_array: PairwiseArrayType, link_cutoff: float) -> List[
	Genotype]:
	"""
		Splits each genotype if any of its members are not related enough to the other members. Genotypes will continue
		to be split until the p-values for all pairwise members are beneath the cutoff.
	Parameters
	----------
	all_genotypes: List[Genotype]
		A list of all genotypes.
	pair_array: PairWiseArrayType
		A mapping of all pairwise p-values.
	link_cutoff: float
		The cuttoff value to determine if a given pair of trajectories is unrelated.

	Returns
	-------
		A list of the new genotypes.

	"""
	for genotype in all_genotypes:
		if len(genotype) > 1:
			combination_pairs = list()
			# Iterate over all possible pairs of genotype members.
			for combination_pair in itertools.combinations(genotype, 2):
				left, right = combination_pair
				# Get the p-value for this pair.
				p_value = get_p_value(left, right, pair_array)
				value = [left, right, p_value.p_value]
				combination_pairs.append(value)
			# Combine all pairs and p-values into a dataframe for convienience.
			genotype_combinations = pandas.DataFrame(combination_pairs, columns = ['left', 'right', 'pvalue'])
			# print(genotype_combinations)
			# Get a dataframe of all trajectories in this genotype which are significantly different than the
			# current pair of trajectories.
			unlinked_trajectories = genotype_combinations[genotype_combinations['pvalue'] <= link_cutoff]

			if len(unlinked_trajectories) != 0:
				# Split the current genotype into two smaller but more internally-related all_genotypes.
				new_genotype_1, new_genotype_2 = split_genotype_in_two(genotype[:], genotype_combinations, pair_array,
					link_cutoff)

				all_genotypes.append(new_genotype_1)
				all_genotypes.append(new_genotype_2)
				all_genotypes.remove(genotype)
	return sorted(all_genotypes, key = lambda s: len(s))


def get_genotypes(timeseries: pandas.DataFrame, options: GenotypeOptions) -> List[Genotype]:
	"""
		Clusters trajectories into genotypes.
	Parameters
	----------
	timeseries: pandas.DataFrame
		The timeseries output from import_timeseries.
		- Columns
			- Population: str
				Name of the population. ex. 'B2'
			- Trajectory: int
				Identifies a unique mutation based on population and posiiton. Should be sorted starting from 1
			- Position: int
				Position of the mutation.
			* timeseries
				The timeseries points will correspond to the timepoints included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
	options: GenotypeOptions
		The options for this script.
	Returns
	-------
	A list of the genotypes, where each genotype is itself a list of the member trajectories.

	"""
	#relative_cutoff = 0.10  # v Calculates the relative similarity between all trajectory pairs.
	#link_cutoff = 0.05  # Calculates the relative similarity between all trajectory pairs.

	# Group all trajectories by the population they belong to. Usually only one population.
	populations = timeseries.groupby(by = 'Population')
	assert populations # Make sure 'populations' is not empty
	all_genotypes = dict()
	for population_id, population_data in populations:
		# Trajectories represent the population frequencies at each timepoint
		# Each row represents a single timepoint, each column represents a mutation.

		# calculate the similarity between all pairs of trajectories in the population.
		pair_array = calculate_pairwise_similarity(
			population_id,
			timeseries,
			detection_cutoff = options.detection_breakpoint,
			fixed_cutoff = options.fixed_breakpoint,
			n_binomial = options.n_binom

		)

		population_genotypes = find_all_genotypes(pair_array, options.similarity_breakpoint)

		# at the end, look at all trajectories that are not listed and
		# append them as their own category.
		# List of all trajectories
		flattened_genotypes = list(itertools.chain.from_iterable(population_genotypes))
		# Any missing trajectories. Will probably be empty
		other_trajectories = [i for i in timeseries['Trajectory'].values if i not in flattened_genotypes]
		population_genotypes.append(other_trajectories)

		# finally, for each genotype, make sure each trajectory pair has some
		# non-trivial linkage (say, >0.0005). this avoids falsely linking together
		# trajectories. if not, divide the offending trajectories into two
		# camps, and sort the rest according to which one they are more
		# closely linked to. repeat until everything is linked.
		# pprint(pair_array)

		while True:
			starting_size_of_the_genotype_array = len(population_genotypes)
			population_genotypes = split_unlinked_genotypes(population_genotypes[:], pair_array, options.difference_breakpoint)
			if len(population_genotypes) == starting_size_of_the_genotype_array:
				break

		all_genotypes[population_id] = population_genotypes
	# TODO Fix so that it returns a genotype for each population
	return population_genotypes


def get_mean_genotypes(all_genotypes: List[Genotype], timeseries: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Calculates the mean frequency of each genotype ate every timepoint.
	Parameters
	----------
	all_genotypes: List[Genotype]
		A list of all genotypes for a given population
	timeseries: pandas.DataFrame
		Must have a 'Trajectory' column along with the columns of the original table that represent timepoints.

	Returns
	-------
	A dataframe where every row corresponds to a genotype (index labels are concatenated member trajectories),
	every column represents a timepoint.
	"""
	mean_genotypes = list()
	for genotype in all_genotypes:
		genotype_timeseries = timeseries[timeseries['Trajectory'].isin(genotype)]
		mean_genotype_timeseries = genotype_timeseries.mean()
		mean_genotype_timeseries.name = "|".join(map(str, genotype))
		mean_genotypes.append(mean_genotype_timeseries)

	mean_genotypes = pandas.DataFrame(mean_genotypes)
	if 'Trajectory' in mean_genotypes:
		mean_genotypes.pop('Trajectory')
	if 'Position' in mean_genotypes:
		mean_genotypes.pop('Position')

	return mean_genotypes


def create_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename'
	)
	parser.add_argument(
		'-o', '--output',
		help = "The folder to save the files to.",
		action = 'store',
		dest = 'output_folder'
	)

	parser.add_argument(
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = 'store',
		dest = 'fixed_breakpoint'
	)
	parser.add_argument(
		"--detected",
		help = "The minimum frequency at which to consider a mutation detected.",
		action = 'store',
		default = 0.03,
		dest = 'detection_breakpoint'
	)

	parser.add_argument(
		"-r", "--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint"
	)

	parser.add_argument(
		"-l", "--difference-cutoff",
		help = "Minimum p-value to consider a pair of genotypes unrelated. Used when splitting genotypes.",
		action = "store",
		default = 0.10,
		dest = "difference_breakpoint"
	)
	parser.add_argument(
		"--matlab",
		help = "Whether to mimic the output of the original matlab script.",
		action = 'store_true',
		dest = 'mode'
	)
	return parser


def workflow(io:Union[Path, pandas.DataFrame], options:GenotypeOptions=None, matlab:bool=True, detection_breakpoint:float=0.03, fixed_breakpoint:float=None)->pandas.DataFrame:
	"""

	Parameters
	----------
	io: Union[Path, pandas.DataFrame]
		Should be either a path to the table of population trajectories or the output of import_timeseries().
	options
	matlab
	detection_breakpoint
	fixed_breakpoint

	Returns
	-------
	pandas.DataFrame
	"""
	if matlab:
		options = GenotypeOptions.from_matlab()
	if options is None:
		if fixed_breakpoint is None:
			fixed_breakpoint = 1-detection_breakpoint

		options = GenotypeOptions.from_breakpoints(detection_breakpoint, fixed_breakpoint)


	if isinstance(io, Path):
		timepoints, info = import_timeseries(io)
	else:
		timepoints = io

	genotypes = get_genotypes(timepoints, options)
	_mean_genotypes = get_mean_genotypes(genotypes, timepoints)

	#_mean_genotypes.to_csv(str(filename.with_suffix('.mean.tsv')), sep = '\t')
	return _mean_genotypes

if __name__ == "__main__":
	cmd_parser = create_parser().parse_args()
	cmd_options = GenotypeOptions.from_parser(cmd_parser)

	workflow(Path(cmd_parser.filename), options = cmd_options)
