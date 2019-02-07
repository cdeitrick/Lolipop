import itertools
from typing import Dict, List, Optional, Tuple

import pandas
from dataclasses import dataclass

try:
	from muller.inheritance.checks import check_additive_background, check_subtractive_background, check_derivative_background
	from muller.options import OrderClusterParameters
except ModuleNotFoundError:
	from .checks import check_additive_background, check_subtractive_background, check_derivative_background
	from options import OrderClusterParameters


@dataclass
class Genotype:
	name: str
	probability: float
	trajectory: pandas.Series
	background: List[str]
	members: List[str]


ClusterType = Dict[str, Genotype]


def add_genotype_background(genotype_label: str, type_genotype: Genotype, nests: Dict[str, Genotype], initial_background_label: str):
	# genotype_label = type_genotype.name
	if genotype_label in nests:
		nests[genotype_label].background += type_genotype.background
		if len(nests[genotype_label].background) > 2 and initial_background_label in nests[genotype_label].background:
			nests[genotype_label].background.remove(initial_background_label)
	else:
		nests[genotype_label] = type_genotype
	return nests


def apply_genotype_checks(type_trajectory: pandas.Series, test_trajectory: pandas.Series, options: OrderClusterParameters) -> Tuple[
	bool, bool, Optional[float]]:
	""" Applies the three checks to `type_trajectory` and `test_trajectory`."""
	additive_check = check_additive_background(
		left = type_trajectory,
		right = test_trajectory,
		double_cutoff = options.additive_background_double_cutoff,
		single_cutoff = options.additive_background_single_cutoff
	)

	subtractive_check = check_subtractive_background(
		left = type_trajectory,
		right = test_trajectory,
		double_cutoff = options.subtractive_background_double_cutoff,
		single_cutoff = options.subtractive_background_single_cutoff
	)

	if subtractive_check:
		delta = None
	else:
		delta = check_derivative_background(
			left = type_trajectory,
			right = test_trajectory,
			detection_cutoff = options.derivative_detection_cutoff
		)

	return additive_check, subtractive_check, delta


def order_clusters(sorted_df: pandas.DataFrame, genotype_members: pandas.Series, options: OrderClusterParameters) -> ClusterType:
	"""
		Orders genotypes by which background they belong to.
	Parameters
	----------
	sorted_df: pandas.DataFrame
		A dataframe of sorted genotypes based on when the genotype was first detected and first fixed.
	genotype_members: pandas.Series
		maps genotypes to their respective members
	options: OrderClusterParameters
		Parameters for the program to use.

	Returns
	-------
	ClusterType
	"""
	# By default the backgrounds should occupy the first n lines of the dataframe
	initial_background = sorted_df.iloc[0]

	nests: Dict[str, Genotype] = {
		initial_background.name: Genotype(
			name = initial_background.name,
			probability = 1.0,
			trajectory = initial_background,
			background = [initial_background.name],
			members = genotype_members.loc[initial_background.name])
	}
	for genotype_label, type_trajectory in sorted_df[1:].iterrows():
		type_members = genotype_members.loc[genotype_label]

		genotype_deltas = list()
		# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
		test_table = sorted_df[:genotype_label].iloc[::-1]
		for test_label, test_trajectory in test_table.iterrows():  # iterate in reverse order
			if genotype_label == test_label:  # The reduced dataframe still contains the genotype being tested.
				continue

			test_background = nests[test_label].background
			type_genotype = Genotype(
				name = genotype_label,
				probability = 1.0,
				trajectory = type_trajectory,
				background = test_background + [genotype_label],
				members = type_members
			)
			additive_check, subtractive_check, delta = apply_genotype_checks(type_trajectory, test_trajectory, options)
			# print(genotype_label, test_label, additive_check, subtractive_check, delta)
			if additive_check:
				nests = add_genotype_background(genotype_label, type_genotype, nests, initial_background.name)
				continue

			if subtractive_check:
				continue

			genotype_deltas.append((test_label, delta))

			if delta > options.derivative_check_cutoff:
				# They are probably on the same background.
				# Need to do one last check: these two genotypes cannot sum to larger than the background.
				nests = add_genotype_background(genotype_label, type_genotype, nests, initial_background.name)
				break

		is_member = genotype_label in itertools.chain.from_iterable(i.background for i in nests.values())
		# Find the test that correlated the most with this genotype
		if genotype_deltas:
			correlated_label, correlated_delta = max(genotype_deltas, key = lambda s: s[1])
		else:
			correlated_label = "N/A"  # Shouldn't end up being used.
			correlated_delta = 0

		if not is_member:
			# if it hasn't been matched with a background, two things can happen
			# (1) IF it is logically possible (i.e., if the sum of it
			# and all existing backgrounds at each time point is ~1 or less),
			# then we can make it its own background. otherwise put
			# it in with the genotype it's most correlated with

			backgrounds = [v.trajectory for k, v in nests.items() if len(v.background) == 1]
			backgrounds = pandas.DataFrame(backgrounds)
			total = backgrounds.sum()
			# Check if There is at least one timepoint where the sum of all current backgrounds is greater than 1
			# or there is at least one timepoint that is less than 1.15.
			exceeds_background_limit = (total > options.new_background_base_cutoff).sum() > 1
			under_cuttoff_limit = (total < options.new_background_significant_cutoff).sum() == 0
			if exceeds_background_limit or not under_cuttoff_limit:
				# add as background
				back = Genotype(
					name = genotype_label,
					probability = 1.0,
					trajectory = type_trajectory,
					background = [genotype_label],
					members = type_members
				)
				nests[genotype_label] = back
			elif correlated_delta > 0:
				correlated_background = nests[correlated_label].background
				type_genotype = Genotype(
					name = genotype_label,
					probability = 0.5,
					trajectory = type_trajectory,
					background = correlated_background + [genotype_label],
					members = type_members
				)
				nests[genotype_label] = type_genotype

			else:
				message = 'SOMETHING HAS GONE HORRIBLY WRONG FOR CLUSTER ' + genotype_label

				raise ValueError(message)

	return nests


if __name__ == "__main__":
	pass
