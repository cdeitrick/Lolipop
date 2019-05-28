import itertools
from typing import Dict, List, Tuple

import pandas

try:
	from muller.clustering.methods.difference import unlink_unrelated_trajectories
	from muller.clustering.metrics.pairwise_calculation_cache import PairwiseCalculationCache
	from muller.clustering.methods.population_genotypes_cache import PopulationGenotypes
except ModuleNotFoundError:
	from .difference import unlink_unrelated_trajectories
	from ..metrics.pairwise_calculation_cache import PairwiseCalculationCache
	from .population_genotypes_cache import PopulationGenotypes


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
	base_genotypes: List[List[str]]
		A set of pre-defined genotypes.

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
	population_genotypes = PopulationGenotypes(genotype_candidates)
	seen = set(itertools.combinations(itertools.chain.from_iterable(genotype_candidates), 2))
	for (left, right), p_value in pairs.items():
		# ignore pairs that have already been sorted into a genotype.
		if (left, right) in seen or (right, left) in seen:
			continue
		seen.add((left, right))
		# are the genotypes related?
		if p_value > relative_cutoff:
			population_genotypes.merge_trajectories(left, right)
	genotype_candidates = population_genotypes.to_list()
	return genotype_candidates


def twostep_method(timeseries: pandas.DataFrame, pair_array: PairwiseCalculationCache, similarity_breakpoint: float, difference_breakpoint: float,
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
	numerical_array = pair_array.asdict()
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

	return [i for i in population_genotypes if i]  # Only return non-empty lists.
