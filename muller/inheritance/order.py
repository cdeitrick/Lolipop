from typing import Dict, List, Tuple

import pandas

try:
	from muller.inheritance.checks import apply_genotype_checks
	from muller.inheritance.cluster import Cluster, Genotype
	from muller.options import OrderClusterParameters
except ModuleNotFoundError:
	from .checks import apply_genotype_checks
	from .cluster import Cluster, Genotype
	from options import OrderClusterParameters


def order_clusters(sorted_df: pandas.DataFrame, genotype_members: pandas.Series, options: OrderClusterParameters) -> Dict[str, Genotype]:
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
	print(sorted_df.to_string())
	initial_background = sorted_df.iloc[0]
	genotype_nests = Cluster(initial_background, timepoints = sorted_df)

	for unnested_label, unnested_trajectory in sorted_df[1:].iterrows():
		genotype_deltas = list()
		# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
		test_table = sorted_df[:unnested_label].iloc[::-1]
		for nested_label, nested_trajectory in test_table.iterrows():  # iterate in reverse order
			if unnested_label == nested_label:  # The reduced dataframe still contains the genotype being tested.
				continue
			test_background = genotype_nests.get(nested_label).background
			type_genotype = Genotype(
				name = unnested_label,
				background = test_background + [unnested_label]
			)
			additive_check, subtractive_check, delta = apply_genotype_checks(unnested_trajectory, nested_trajectory, options)
			# print(genotype_label, test_label, additive_check, subtractive_check, delta)
			if additive_check:
				genotype_nests.add_genotype_to_background(unnested_label, type_genotype)
				continue

			if subtractive_check:
				continue

			genotype_deltas.append((nested_label, delta))

			if delta > options.derivative_check_cutoff:
				# They are probably on the same background.
				# Need to do one last check: these two genotypes cannot sum to larger than the background.
				genotype_nests.add_genotype_to_background(unnested_label, type_genotype)
				break

		is_member = genotype_nests.is_a_member(unnested_label)
		if not is_member:
			print(unnested_label)
			# if it hasn't been matched with a background, two things can happen
			# (1) IF it is logically possible (i.e., if the sum of it
			# and all existing backgrounds at each time point is ~1 or less),
			# then we can make it its own background. otherwise put
			# it in with the genotype it's most correlated with
			# Find the test that correlated the most with this genotype
			result = background_heuristic(genotype_nests, genotype_deltas, unnested_trajectory, options)
			if result:
				genotype_nests.add_genotype_to_background(unnested_label, result)
	return genotype_nests.nests


def get_maximum_genotype_delta(genotype_deltas: List[Tuple[str, float]]) -> Tuple[str, float]:
	if genotype_deltas:
		correlated_label, correlated_delta = max(genotype_deltas, key = lambda s: s[1])
	else:
		correlated_label = "N/A"  # Shouldn't end up being used.
		correlated_delta = 0
	return correlated_label, correlated_delta


# noinspection PyUnresolvedReferences,PyTypeChecker
def background_heuristic(genotype_nests: Cluster, genotype_deltas: List[Tuple[str, float]], unnested_trajectory: pandas.Series,
		options: OrderClusterParameters):
	unnested_label = unnested_trajectory.name
	back = Genotype(
		name = unnested_label,
		background = [unnested_label]
	)

	correlated_label, correlated_delta = get_maximum_genotype_delta(genotype_deltas)
	current_background_total = genotype_nests.get_sum_of_backgrounds()

	try:
		correlated_background = genotype_nests.get(correlated_label).background
	except KeyError:
		correlated_background = None

	exceeds_background_limit = (current_background_total > options.new_background_base_cutoff).sum() > 1
	under_cuttoff_limit = (current_background_total < options.new_background_significant_cutoff).sum() == 0

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
	from pprint import pprint
	from import_data import import_table_from_string
	string = """
		Genotype	0	17	25	44	66	75	90
		genotype-1	0	0	0.261	1	1	1	1
		genotype-6	0	0	0	0.273	0.781	1	1
		genotype-3	0	0	0	0	0	1	1
		genotype-4	0	0	0	0.525	0.454	0.911	0.91
		genotype-5	0	0	0	0.147	0.45	0.924	0.887
		genotype-11	0	0	0	0	0.278	0.822	0.803
		genotype-2	0	0.38	0.432	0	0	0	0
		genotype-8	0	0	0	0.403	0.489	0.057	0.08
		genotype-14	0	0	0	0	0	0.2675	0.326
		genotype-10	0	0	0	0.138	0.295	0	0.081
		genotype-12	0	0	0	0	0.2335	0.133	0.0375
		genotype-7	0	0	0	0.188	0.171	0.232	0.244
		genotype-9	0	0	0.117	0	0	0	0.103
		genotype-13	0	0	0.033	0.106	0.1065	0	0
		genotype-15	0	0	0	0.1145	0	0.1205	0.0615
	"""
	table = import_table_from_string(string, index = 'Genotype')
	result = order_clusters(table, None, OrderClusterParameters.from_breakpoints(.03, .97))
	pprint(result)
