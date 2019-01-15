from typing import List, Tuple, Dict
import itertools
import pandas

def _divide_genotype(genotype: List[str], unlinked_trajectories: pandas.DataFrame,
		pair_array: Dict[Tuple[str,str], float], link_cut: float) -> Tuple[List[str], List[str]]:
	"""
		Splits a genotype into smaller muller_genotypes if some members are not related to some other members. This may happen
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
		a tuple of two muller_genotypes.
	"""
	# Find the index of the minimum p-value after subtracting the link cutoff.
	minimum_pvalue_index = (unlinked_trajectories['pvalue'] - link_cut).abs().idxmin()
	# Get the row with the identified minimu mp-value.
	minimum_pvalue_genotype = unlinked_trajectories.loc[minimum_pvalue_index]

	# Form two new muller_genotypes based on the two trajectories corresponding to the minimum p-value
	new_genotype_1_base = minimum_pvalue_genotype['left']
	new_genotype_2_base = minimum_pvalue_genotype['right']

	# Make sure genotype 1 includes the lower id-value. Not important, but maintains parity with matlab script.
	new_genotype_1_base, new_genotype_2_base = sorted([new_genotype_1_base, new_genotype_2_base])
	new_genotype_1 = [new_genotype_1_base]
	new_genotype_2 = [new_genotype_2_base]

	# Sort each member in the current genotype into one of the new muller_genotypes.
	for genotype_member in genotype:
		genotype_member = genotype_member
		# Check if the current genotype member is already contained in one of the muller_genotypes.
		# Should only be one of the two trajectories used to form a new genotype.
		if genotype_member in new_genotype_1 or genotype_member in new_genotype_2:
			pass
		else:
			# Use the highest p-value to determine which genotype to add the member to.
			# P-values should correspond to the current member and the base member of the new muller_genotypes.
			p_value_1 = pair_array[new_genotype_1_base, genotype_member]
			p_value_2 = pair_array[new_genotype_2_base, genotype_member]

			if p_value_1 >= p_value_2:
				new_genotype_1.append(genotype_member)
			else:
				new_genotype_2.append(genotype_member)
	return new_genotype_1, new_genotype_2

def unlink_unrelated_trajectories(all_genotypes: List[List[str]], pair_array: Dict[Tuple[str,str], float], link_cutoff: float) -> List[List[str]]:
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
			combination_pairs = [(left, right, pair_array[left, right]) for left, right in itertools.combinations(genotype, 2)]
			# Combine all pairs and p-values into a dataframe for convienience.
			genotype_combinations = pandas.DataFrame(combination_pairs, columns = ['left', 'right', 'pvalue'])
			# Get a dataframe of all trajectories in this genotype which are significantly different than the
			# current pair of trajectories.
			unlinked_trajectories = genotype_combinations[genotype_combinations['pvalue'] < link_cutoff]

			if len(unlinked_trajectories) != 0:
				# Split the current genotype into two smaller but more internally-related all_genotypes.
				new_genotype_1, new_genotype_2 = _divide_genotype(genotype[:], genotype_combinations, pair_array, link_cutoff)

				all_genotypes.append(new_genotype_1)
				all_genotypes.append(new_genotype_2)
				all_genotypes.remove(genotype)
	# noinspection PyTypeChecker
	return sorted(all_genotypes, key = len)
