import itertools
from typing import Dict, List, Tuple

import pandas
from dataclasses import dataclass


@dataclass
class Genotype:
	name: str
	probability: float
	trajectory: pandas.Series
	background: List[str]
	members: List[str]


ClusterType = Dict[str, Genotype]


@dataclass
class OrderClusterParameters:
	additive_background_double_cutoff: float
	additive_background_single_cutoff: float
	subtractive_background_double_cutoff: float
	subtractive_background_single_cutoff: float
	derivative_detection_cutoff: float
	derivative_check_cutoff: float
	new_background_base_cutoff: float
	new_background_significant_cutoff: float

	@classmethod
	def from_matlab(cls) -> 'OrderClusterParameters':
		return OrderClusterParameters(
			additive_background_double_cutoff = 1.03,
			additive_background_single_cutoff = 1.15,
			subtractive_background_double_cutoff = -0.02,
			subtractive_background_single_cutoff = -0.15,
			derivative_detection_cutoff = 0.02,
			derivative_check_cutoff = 0.01,
			new_background_base_cutoff = 1.0,
			new_background_significant_cutoff = 1.15
		)

	@classmethod
	def from_breakpoints(cls, detection_breakpoint: float, significant_breakpoint: float) -> 'OrderClusterParameters':
		return OrderClusterParameters(
			additive_background_double_cutoff = 1 + detection_breakpoint,
			additive_background_single_cutoff = 1 + significant_breakpoint,
			subtractive_background_double_cutoff = -detection_breakpoint,
			subtractive_background_single_cutoff = -significant_breakpoint,
			derivative_check_cutoff = 0.01,
			derivative_detection_cutoff = detection_breakpoint,
			new_background_base_cutoff = 1 + detection_breakpoint,
			new_background_significant_cutoff = 1 + significant_breakpoint
		)


def check_additive_background(left: pandas.Series, right: pandas.Series, double_cutoff: float,
		single_cutoff: float) -> bool:
	# Want to check if a genotype's sum of frequencies is consistently greater than 1.0
	trajectorysum = right + left
	double_check = (trajectorysum > double_cutoff).sum() > 0  # Implicit conversion from bool to int.
	single_check = (trajectorysum > single_cutoff).sum() > 0
	return double_check or single_check


def check_subtractive_background(left: pandas.Series, right: pandas.Series, double_cutoff: float,
		single_cutoff: float) -> bool:
	# Check if the current genotype is over 15% larger than the other.
	diff_trajectory = right - left
	double_diff_trajectory = (diff_trajectory < double_cutoff).sum() > 1
	single_diff_trajectory = (diff_trajectory < single_cutoff).sum() > 0  # implicit conversion from bool to int
	return double_diff_trajectory or single_diff_trajectory


# noinspection PyTypeChecker
def check_derivative_background(left: pandas.Series, right: pandas.Series, detection_cutoff: float) -> float:
	# Look at the first point when both trajectories are non-zero, and observe delta over time.
	# Are derivatives correlated or anti-correlated? If sufficiently correlated, they are on the same background.
	right_trajectory_filtered: pandas.Series = right > detection_cutoff  # All positions > 0.02 will evaluate to True
	left_trajectory_filtered: pandas.Series = left > detection_cutoff
	# Find the points where neither left nor right is zero.
	startpoint = max([right_trajectory_filtered.idxmax(), left_trajectory_filtered.idxmax()])
	endpoint = min([right_trajectory_filtered[::-1].idxmax(), left_trajectory_filtered[::-1].idxmax()])

	# Convert the index labels to their corresponding integer position in the index.
	# Pandas can't use labels to get a slice of a series.
	start_index = left.index.get_loc(startpoint)
	end_index = left.index.get_loc(endpoint)
	# The matlab script uses the difference of the following timepoint versus the curent one,
	# while pandas default behaviour is current timepoint - previous timepoint.
	# Need to indicate difference of following timepoint, and reverse polarity(negative sign)
	right_derivative = right.diff(-1)
	left_derivative = left.diff(-1)

	right_window = right_derivative[start_index:end_index] * -1
	left_window = left_derivative[start_index:end_index] * -1
	delta = left_window.dot(right_window)

	return delta


def add_genotype_bakground(genotype_label: str, type_genotype: Genotype, nests: Dict[str, Genotype], initial_background_label: str):
	# genotype_label = type_genotype.name
	if genotype_label in nests:
		nests[genotype_label].background += type_genotype.background
		if len(nests[genotype_label].background) > 2 and initial_background_label in nests[genotype_label].background:
			nests[genotype_label].background.remove(initial_background_label)
	else:
		nests[genotype_label] = type_genotype
	return nests


def apply_genotype_checks(type_trajectory: pandas.Series, test_trajectory: pandas.Series, options: OrderClusterParameters) -> Tuple[
	bool, bool, float]:
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
				nests = add_genotype_bakground(genotype_label, type_genotype, nests, initial_background.name)
				continue

			if subtractive_check:
				continue

			genotype_deltas.append((test_label, delta))

			if delta > options.derivative_check_cutoff:
				# They are probably on the same background.
				# Need to do one last check: these two genotypes cannot sum to larger than the background.
				nests = add_genotype_bakground(genotype_label, type_genotype, nests, initial_background.name)
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
	index = [0, 1, 3, 4, 6, 7, 9, 10, 12]
	left = pandas.Series([0, 0, 0.47, 1, 1, 1, 1, 0.173, 0.169], index = index)
	right = pandas.Series([0, 0, 0, 0, 0, 0.347, 0.449, 0, 0.097], index = index)
	additive_check = check_additive_background(left, right, .03, .15)
	subtractive_check = check_subtractive_background(left, right, -.03, -.15)
	derivative_check = check_derivative_background(left, right, .03)

	print(additive_check, subtractive_check, derivative_check)
