import itertools
from typing import Dict, List, Tuple

import pandas


# noinspection PyUnresolvedReferences,PyTypeChecker
def _find_weakest_pair(trajectories: pandas.Series, link_cut: float) -> Tuple[str, str]:
	return (trajectories - link_cut).abs().idxmin()


# noinspection PyTypeChecker
def _divide_genotype(genotype: List[str], unlinked_trajectories: pandas.Series, link_cut: float) -> Tuple[List[str], List[str]]:
	"""
		Splits a genotype into smaller muller_genotypes if some members are not related to some other members. This may happen
		when a member is included into a genotype due to its paired member but ends up not related to some of
		the other members that where found later.
	Parameters
	----------
	genotype: Genotype
		The genotype to split.
	unlinked_trajectories: pandas.Series
		A series of trajectory pairs and their corresponding p-value.
		Index:
			- left
			- right
		Values:
			- pvalue
	link_cut: float
		The cuttoff value to choose whether a member is related to the other members.

	Returns
	-------
		a tuple of two muller_genotypes.
	"""
	# Find the index of the minimum p-value after subtracting the link cutoff.
	# Form two new muller_genotypes based on the two trajectories corresponding to the minimum p-value
	# noinspection PyUnresolvedReferences
	# new_genotype_1_base, new_genotype_2_base = (unlinked_trajectories - link_cut).abs().idxmin()
	new_genotype_1_base, new_genotype_2_base = _find_weakest_pair(unlinked_trajectories, link_cut)
	# Get the row with the identified minimum mp-value.

	# Make sure genotype 1 includes the lower id-value. Not important, but maintains parity with matlab script.
	new_genotype_1_base, new_genotype_2_base = sorted([new_genotype_1_base, new_genotype_2_base])
	new_genotype_1 = [new_genotype_1_base]
	new_genotype_2 = [new_genotype_2_base]

	# Sort each member in the current genotype into one of the new muller_genotypes.
	for genotype_member in genotype:
		# Check if the current genotype member is already contained in one of the muller_genotypes.
		# Should only be one of the two trajectories used to form a new genotype.
		if not (genotype_member in new_genotype_1 or genotype_member in new_genotype_2):
			# Use the highest p-value to determine which genotype to add the member to.
			# P-values should correspond to the current member and the base member of the new muller_genotypes.
			p_value_1 = unlinked_trajectories.loc[new_genotype_1_base, genotype_member]
			p_value_2 = unlinked_trajectories.loc[new_genotype_2_base, genotype_member]

			if p_value_1 >= p_value_2:
				new_genotype_1.append(genotype_member)
			else:
				new_genotype_2.append(genotype_member)
	return new_genotype_1, new_genotype_2


def unlink_unrelated_trajectories(all_genotypes: List[List[str]], pair_array: Dict[Tuple[str, str], float], link_cutoff: float) -> List[List[str]]:
	"""
		Splits each genotype if any of its members are not related enough to the other members. Genotypes will continue
		to be split until the p-values for all pairwise members are beneath the cutoff.
	Parameters
	----------
	all_genotypes: List[Genotype]
		A list of all muller_genotypes.
	pair_array: PairWiseArrayType
		A mapping of all pairwise p-values.
	link_cutoff: float
		The cuttoff value to determine if a given pair of trajectories is unrelated.

	Returns
	-------
		A list of the new muller_genotypes.

	"""
	for genotype in all_genotypes:
		if len(genotype) > 1:
			# Iterate over all possible pairs of genotype members.
			combination_pairs = [(left, right, pair_array[left, right]) for left, right in itertools.permutations(genotype, 2)]
			# Combine all pairs and p-values into a dataframe for convienience.
			genotype_combinations: pandas.DataFrame = pandas.DataFrame(combination_pairs, columns = ['left', 'right', 'pvalue'])
			genotype_combinations: pandas.Series = genotype_combinations.set_index(['left', 'right'])['pvalue']
			# Get a dataframe of all trajectories in this genotype which are significantly different than the
			# current pair of trajectories.
			# noinspection PyTypeChecker
			unlinked_trajectories = genotype_combinations[genotype_combinations < link_cutoff].dropna()

			if len(unlinked_trajectories) != 0:
				# Split the current genotype into two smaller but more internally-related all_genotypes.
				new_genotype_1, new_genotype_2 = _divide_genotype(genotype[:], genotype_combinations, link_cutoff)

				all_genotypes.append(new_genotype_1)
				all_genotypes.append(new_genotype_2)
				all_genotypes.remove(genotype)
	# noinspection PyTypeChecker
	return sorted(all_genotypes, key = len)
