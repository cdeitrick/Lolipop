import argparse
import itertools
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas
from dataclasses import dataclass

try:
	from muller.import_table import import_trajectory_table
except ModuleNotFoundError:
	# noinspection PyUnresolvedReferences
	from import_table import import_trajectory_table


@dataclass
class PairArrayValue:
	"""
		Holds the values assigned to parray in the original script.
	"""
	left: str
	right: str
	p_value: float
	sigma_value: float
	difbar: float


# PairwiseArrayType = Dict[Tuple[str, str], PairArrayValue]
PairwiseArrayType = Dict[Tuple[str, str], float]
PAIRWISE_P_VALUES: PairwiseArrayType = None
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


# noinspection PyTypeChecker
def calculate_p_value(left: pandas.Series, right: pandas.Series, detected_cutoff: float, fixed_cutoff: float) -> float:
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
	left, right: pandas.Series
		- index: int
			Timepoints
		- values: float
			Frequencies
	detected_cutoff: float
	fixed_cutoff: float


	Returns
	-------
	float
	"""
	# Merge into a dataframe for convienience
	df = pandas.concat([left, right], axis = 1)

	# Remove timepoints where at least one trajectory was not fixed or undetected.
	not_fixed_df = df[(df < fixed_cutoff).any(axis = 1)]

	not_detected_fixed_df = not_fixed_df[(not_fixed_df > detected_cutoff).any(axis = 1)]

	if not_detected_fixed_df.empty:
		left_fixed: pandas.Series = left[left > fixed_cutoff]
		right_fixed: pandas.Series = right[right > fixed_cutoff]

		if left_fixed.empty and right_fixed.empty:
			# Both are undetected
			p_value = 1
		else:
			overlap = set(left_fixed.index) and set(right_fixed.index)
			p_value = int(len(left_fixed) > 2 and len(right_fixed) > 2 and len(overlap) > 2)

	else:
		# Find the mean frequency of each timepoint
		# index is timepoints,  values are frequencies
		mean: pandas.Series = not_detected_fixed_df.mean(axis = 1)

		n_binom = len(mean)  # WARNING: is not compatible with matlab scripts for n != 5
		# Calculate sigma_freq
		sigma_freq: pandas.Series = (mean * (1 - mean)) / n_binom
		# Difference of frequencies at each timepoint
		difference: pandas.Series = not_detected_fixed_df.iloc[:, 0] - not_detected_fixed_df.iloc[:, 1]
		sigma_pair: float = math.sqrt(sigma_freq.sum()) / len(difference)
		# Sum of differences
		difference_mean: float = abs(difference).sum() / len(difference)

		X = difference_mean / (math.sqrt(2) * sigma_pair)

		p_value: float = 1 - math.erf(X)

	return p_value


def calculate_trajectory_similarity(trajectories: pandas.DataFrame, detection_cutoff: float,
		fixed_cutoff: float) -> PairwiseArrayType:
	"""
	Parameters
	----------
	trajectories: pandas.DataFrame
		A table of mutational trajectories.
	detection_cutoff: float
	fixed_cutoff: float

	Returns
	-------
	dict of PairArrayValue
	Each key in the dictionary corresponds to a pair of trajectory ids which map to the p-value for that pair.
	The order of ids does not matter.
	"""
	combos: List[Tuple[str, str]] = sorted(itertools.combinations(trajectories.index, 2))
	pair_array = dict()
	for left, right in combos:
		left_trajectories = trajectories.loc[left]
		right_trajectories = trajectories.loc[right]

		p_value = calculate_p_value(left_trajectories, right_trajectories, detection_cutoff, fixed_cutoff)

		pair_array[(left, right)] = p_value
		pair_array[(right, left)] = p_value

	return pair_array


def _find_genotype_from_trajectory(element: str, all_genotypes: List[List[str]]) -> Optional[List[str]]:
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


def find_all_genotypes(pairs: PairwiseArrayType, relative_cutoff: float) -> List[List[str]]:
	"""
		Clusters all trajectories into related genotypes.
		By default the first trajectory makes a genotype category
	Parameters
	----------
	pairs: PairwiseArrayType
		A dictionary mapping pairs to p-values.
	relative_cutoff: float; default 0.05
		The cutoff indicating two trajectories are related.

	Returns
	-------
		A list of all genotypes (each of which is a list of ints)
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
					# Remove the redundant genotypes.
					genotype_candidates.remove(genotype_right)

			elif genotype_left:
				genotype_left.append(right)
			elif genotype_right:
				genotype_right.append(left)
			else:
				# Neither element is listed. Create a new genotype
				genotype_candidates.append([left, right])
	return genotype_candidates


def _divide_genotype(genotype: List[str], unlinked_trajectories: pandas.DataFrame,
		pair_array: PairwiseArrayType, link_cut: float) -> Tuple[List[str], List[str]]:
	"""
		Splits a genotype into smaller genotypes if some members are not related to some other members. This may happen
		when a member is included into a genotype due to its paired member but ends up not related to some of
		the other members that where found later.
	Parameters
	----------
	genotype: Genotype
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
	new_genotype_1_base = minimum_pvalue_genotype['left']
	new_genotype_2_base = minimum_pvalue_genotype['right']

	# Make sure genotype 1 includes the lower id-value. Not important, but maintains parity with matlab script.
	new_genotype_1_base, new_genotype_2_base = sorted([new_genotype_1_base, new_genotype_2_base])
	new_genotype_1 = [new_genotype_1_base]
	new_genotype_2 = [new_genotype_2_base]

	# Sort each member in the current genotype into one of the new genotypes.
	for genotype_member in genotype:
		genotype_member = genotype_member
		# Check if the current genotype member is already contained in one of the genotypes.
		# Should only be one of the two trajectories used to form a new genotype.
		if genotype_member in new_genotype_1 or genotype_member in new_genotype_2:
			pass
		else:
			# Use the highest p-value to determine which genotype to add the member to.
			# P-values should correspond to the current member and the base member of the new genotypes.
			p_value_1 = pair_array[new_genotype_1_base, genotype_member]
			p_value_2 = pair_array[new_genotype_2_base, genotype_member]

			if p_value_1 >= p_value_2:
				new_genotype_1.append(genotype_member)
			else:
				new_genotype_2.append(genotype_member)
	return new_genotype_1, new_genotype_2


def _unlink_unrelated_trajectories(all_genotypes: List[List[str]], pair_array: PairwiseArrayType, link_cutoff: float) -> List[
	List[str]]:
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
				p_value = pair_array[left, right]
				value = [left, right, p_value]
				combination_pairs.append(value)
			# Combine all pairs and p-values into a dataframe for convienience.
			genotype_combinations = pandas.DataFrame(combination_pairs, columns = ['left', 'right', 'pvalue'])
			# print(genotype_combinations)
			# Get a dataframe of all trajectories in this genotype which are significantly different than the
			# current pair of trajectories.
			unlinked_trajectories = genotype_combinations[genotype_combinations['pvalue'] <= link_cutoff]

			if len(unlinked_trajectories) != 0:
				# Split the current genotype into two smaller but more internally-related all_genotypes.
				new_genotype_1, new_genotype_2 = _divide_genotype(genotype[:], genotype_combinations, pair_array,
					link_cutoff)

				all_genotypes.append(new_genotype_1)
				all_genotypes.append(new_genotype_2)
				all_genotypes.remove(genotype)
	return sorted(all_genotypes, key = len)


def calculate_genotypes_from_population(timeseries: pandas.DataFrame, options: GenotypeOptions) -> List[List[str]]:
	"""
		Clusters trajectories into genotypes.
	Parameters
	----------
	timeseries: pandas.DataFrame
		The timeseries output from import_timeseries.
		- Columns
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

	# Trajectories represent the population frequencies at each timepoint
	# Each row represents a single timepoint, each column represents a mutation.

	# calculate the similarity between all pairs of trajectories in the population.
	global PAIRWISE_P_VALUES

	if PAIRWISE_P_VALUES:
		for key in sorted(PAIRWISE_P_VALUES.keys()):
			left_key, right_key = key
			if left_key not in timeseries.index or right_key not in timeseries.index:
				PAIRWISE_P_VALUES.pop(key)
		pair_array = PAIRWISE_P_VALUES

	else:
		pair_array = calculate_trajectory_similarity(
			timeseries,
			detection_cutoff = options.detection_breakpoint,
			fixed_cutoff = options.fixed_breakpoint
		)
		PAIRWISE_P_VALUES = pair_array

	population_genotypes = find_all_genotypes(pair_array, options.similarity_breakpoint)

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
		population_genotypes = _unlink_unrelated_trajectories(population_genotypes[:], pair_array,
			options.difference_breakpoint)
		if len(population_genotypes) == starting_size_of_the_genotype_array:
			break

	# all_genotypes[population_id] = population_genotypes
	# TODO Fix so that it returns a genotype for each population
	return [i for i in population_genotypes if i]  # Only return non-empty lists.


def calculate_mean_genotype(all_genotypes: List[List[str]], timeseries: pandas.DataFrame) -> pandas.DataFrame:
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
	A dataframe where every row corresponds to a genotype.
	member trajectories are listed under the 'members' column.
	every column represents a timepoint.
	"""

	mean_genotypes = list()

	for index, genotype in enumerate(all_genotypes, start = 1):
		genotype_timeseries = timeseries.loc[genotype]
		mean_genotype_timeseries = genotype_timeseries.mean()
		mean_genotype_timeseries['members'] = "|".join(map(str, genotype))
		mean_genotype_timeseries.name = "genotype-{}".format(index)
		mean_genotypes.append(mean_genotype_timeseries)

	mean_genotypes = pandas.DataFrame(mean_genotypes)

	if 'Trajectory' in mean_genotypes:
		mean_genotypes.pop('Trajectory')
	if 'Position' in mean_genotypes:
		mean_genotypes.pop('Position')

	if 'members' in mean_genotypes.columns:
		members = mean_genotypes.pop('members')
	else:
		members = None

	mean_genotypes = mean_genotypes[sorted(mean_genotypes.columns)]

	if members is not None:
		mean_genotypes['members'] = members
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


def workflow(io: Union[Path, pandas.DataFrame], options: GenotypeOptions = None, matlab: bool = False,
		detection_breakpoint: float = 0.03, fixed_breakpoint: float = None) -> pandas.DataFrame:
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
			fixed_breakpoint = 1 - detection_breakpoint

		options = GenotypeOptions.from_breakpoints(detection_breakpoint, fixed_breakpoint)

	if isinstance(io, Path):
		timepoints, _ = import_trajectory_table(io)
	else:
		timepoints = io

	genotypes = calculate_genotypes_from_population(timepoints, options)

	_mean_genotypes = calculate_mean_genotype(genotypes, timepoints)

	# _mean_genotypes.to_csv(str(filename.with_suffix('.mean.tsv')), sep = '\t')
	return _mean_genotypes


if __name__ == "__main__":
	pass