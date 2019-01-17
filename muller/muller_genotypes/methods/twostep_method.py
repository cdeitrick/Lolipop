import itertools
from typing import Dict, List, Optional, Tuple

import pandas

try:
	from muller_genotypes.methods.difference import unlink_unrelated_trajectories
	from muller_genotypes.metrics.pairwise_calculation import PairwiseCalculation
except ModuleNotFoundError:
	from .difference import unlink_unrelated_trajectories
	from ..metrics.pairwise_calculation import PairwiseCalculation


def _group_trajectories_into_genotypes(pairs: Dict[Tuple[str, str], float], relative_cutoff: float, base_genotypes: List[List[str]] = None) -> List[
	List[str]]:
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
	if base_genotypes:
		genotype_candidates = base_genotypes
	elif pairs:
		genotype_candidates = [[min(pairs.keys())[0]]]  # by default the first trajectory forms the first genotype.
	else:
		genotype_candidates = []

	seen = set()
	for (left, right), p_value in pairs.items():
		# ignore pairs that have already been sorted into a genotype.
		if (left, right) in seen or (right, left) in seen:
			continue
		seen.add((left, right))
		# are the genotypes related?
		if p_value > relative_cutoff:
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


def _find_genotype_from_trajectory(element: str, all_genotypes: List[List[str]]) -> Optional[List[str]]:
	"""
		Finds the genotype that contains the trajectory.
	Parameters
	----------
	element: str
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


def twostep_method(timeseries: pandas.DataFrame, pair_array: PairwiseCalculation, similarity_breakpoint: float, difference_breakpoint: float,
		starting_genotypes: List[List[str]]) -> List[List[str]]:
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
	pair_array: PairwiseArrayType
	similarity_breakpoint: float
	difference_breakpoint: float
	starting_genotypes: List[List[str]]
		A list of trajectories to use as the initial genotypes. Each sub-list represents a single genotype.
	Returns
	-------
	List[List[str]]
		A list of the genotypes, where each genotype is itself a list of the name of each member trajectory.
		Ex. [
			[A1, A2, A3],
			[B1, B2],
			[C1]
		]
	"""

	# Trajectories represent the population frequencies at each timepoint
	# Each row represents a single timepoint, each column represents a mutation.
	numerical_array = pair_array.asitem('pvalue')
	population_genotypes = _group_trajectories_into_genotypes(numerical_array, similarity_breakpoint, starting_genotypes)

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
		population_genotypes = unlink_unrelated_trajectories(population_genotypes[:], numerical_array, difference_breakpoint)
		if len(population_genotypes) == starting_size_of_the_genotype_array:
			break

	# all_genotypes[population_id] = population_genotypes
	# TODO Fix so that it returns a genotype for each population
	return [i for i in population_genotypes if i]  # Only return non-empty lists.
