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

class Cluster:
	def __init__(self, initial_background:pandas.Series, initial_members:str):
		self.initial_background_label:str = initial_background.name

		initial_genotype = Genotype(
			name = initial_background.name,
			probability = 1.0,
			trajectory = initial_background,
			background = [initial_background.name],
			members = [initial_members]
		)

		self.nests = {initial_background.name:initial_genotype}

	def add_genotype_to_background(self, genotype_label:str, type_genotype:Genotype)->None:
		if genotype_label in self.nests:
			self.nests[genotype_label].background += type_genotype.background

			if len(self.nests[genotype_label].background) > 2 and self.initial_background_label in self.nests[genotype_label].background:
				self.nests[genotype_label].background.remove(self.initial_background_label)

		else:
			self.nests[genotype_label] = type_genotype



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
	initial_members = genotype_members.loc[initial_background.name]
	nests: Dict[str, Genotype] = {
		initial_background.name: Genotype(
			name = initial_background.name,
			probability = 1.0,
			trajectory = initial_background,
			background = [initial_background.name],
			members = genotype_members.loc[initial_background.name])
	}
	genotype_nests: Cluster(initial_background, initial_members)

	for unnested_label, unnested_trajectory in sorted_df[1:].iterrows():
		type_members = genotype_members.loc[unnested_label]
		genotype_nests[unnested_label] = []
		genotype_deltas = list()
		# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
		test_table = sorted_df[:unnested_label].iloc[::-1]
		for nested_label, nested_trajectory in test_table.iterrows():  # iterate in reverse order
			if unnested_label == nested_label:  # The reduced dataframe still contains the genotype being tested.
				continue

			test_background = nests[nested_label].background
			type_genotype = Genotype(
				name = unnested_label,
				probability = 1.0,
				trajectory = unnested_trajectory,
				background = test_background + [unnested_label],
				members = type_members
			)
			additive_check, subtractive_check, delta = apply_genotype_checks(unnested_trajectory, nested_trajectory, options)
			# print(genotype_label, test_label, additive_check, subtractive_check, delta)
			if additive_check:
				nests = add_genotype_background(unnested_label, type_genotype, nests, initial_background.name)
				genotype_nests[unnested_label].append(nested_label)
				print("Add ", unnested_label, nested_label)
				continue

			if subtractive_check:
				continue

			genotype_deltas.append((nested_label, delta))

			if delta > options.derivative_check_cutoff:
				# They are probably on the same background.
				# Need to do one last check: these two genotypes cannot sum to larger than the background.
				nests = add_genotype_background(unnested_label, type_genotype, nests, initial_background.name)
				genotype_nests[unnested_label].append(nested_label)
				print("del ", unnested_label, nested_label)
				break

		is_member = unnested_label in itertools.chain.from_iterable(i.background for i in nests.values())

		if not is_member:
			# if it hasn't been matched with a background, two things can happen
			# (1) IF it is logically possible (i.e., if the sum of it
			# and all existing backgrounds at each time point is ~1 or less),
			# then we can make it its own background. otherwise put
			# it in with the genotype it's most correlated with
			# Find the test that correlated the most with this genotype
			back = Genotype(
				name = unnested_label,
				probability = 1.0,
				trajectory = unnested_trajectory,
				background = [unnested_label],
				members = type_members
			)
			correlated_label, correlated_delta = get_maximum_genotype_delta(genotype_deltas)
			current_background_total = get_sum_of_backgrounds(nests)
			correlated_background = nests[correlated_label].background
			result = background_heuristic(
				back,
				current_background_total,
				correlated_delta,
				correlated_background,
				options.new_background_base_cutoff,
				options.new_background_significant_cutoff
			)
			if result:
				nests[unnested_label] = result
	from pprint import pprint
	pprint(genotype_nests, width = 300)
	return nests


def get_maximum_genotype_delta(genotype_deltas: List[Tuple[str, float]]) -> Tuple[str, float]:
	if genotype_deltas:
		correlated_label, correlated_delta = max(genotype_deltas, key = lambda s: s[1])
	else:
		correlated_label = "N/A"  # Shouldn't end up being used.
		correlated_delta = 0
	return correlated_label, correlated_delta


def is_a_background(element: Genotype) -> bool:
	return len(element.background) == 1


def get_sum_of_backgrounds(nests: Dict[str, Genotype]) -> pandas.Series:
	backgrounds = [v.trajectory for k, v in nests.items() if is_a_background(v)]
	backgrounds = pandas.DataFrame(backgrounds)
	total = backgrounds.sum()
	return total


def background_heuristic(back: Genotype, background_total: pandas.Series, correlated_delta, correlated_background,
		base_cutoff: float, sig_cutoff: float):
	# Check if There is at least one timepoint where the sum of all current backgrounds is greater than 1
	# or there is at least one timepoint that is less than 1.15.
	exceeds_background_limit = (background_total > base_cutoff).sum() > 1
	under_cuttoff_limit = (background_total < sig_cutoff).sum() == 0

	if exceeds_background_limit or not under_cuttoff_limit:
		# add as background
		result = back
	elif correlated_delta > 0:
		back.background = correlated_background + [back.background]
		result = back
	else:
		# message = 'SOMETHING HAS GONE HORRIBLY WRONG FOR CLUSTER ' + unnested_label
		result = None
	return result


if __name__ == "__main__":
	from import_data import import_table_from_string

	sorted_string = """
		Genotype	0	17	25	44	66	75	90
		1	0	0	0.261	1	1	1	1
		7	0	0	0	0.273	0.781	1	1
		6	0	0	0	0	0	1	1
		2	0	0	0	0.525	0.454	0.911	0.91
		3	0	0	0	0.147	0.45	0.924	0.887
		14	0	0.38	0.432	0	0	0	0
		9	0	0	0	0	0	0.269	0.34
		17	0	0	0	0	0	0.266	0.312
		20	0	0	0	0.138	0.295	0	0.081
		13	0	0	0	0	0.258	0.057	0.075
		16	0	0	0	0	0.209	0.209	0
		10	0	0	0.117	0	0	0	0.103
		15	0	0	0.066	0.104	0.062	0	0
		11	0	0	0	0.108	0.151	0	0
	"""
	options = OrderClusterParameters.from_breakpoints(0.03, 0.97)
	table = import_table_from_string(sorted_string, index = 'Genotype')
	table['members'] = ""
	members = table.pop('members')

	order_clusters(table, members, options)
