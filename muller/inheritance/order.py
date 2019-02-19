import logging
from typing import Dict, List, Tuple

import pandas

logger = logging.getLogger(__file__)
try:
	from muller.inheritance.checks import apply_genotype_checks
	from muller.inheritance import checks
	from muller.inheritance.cluster import Cluster
	from muller.options import OrderClusterParameters
except ModuleNotFoundError:
	from . import checks
	from .checks import apply_genotype_checks
	from .cluster import Cluster
	from options import OrderClusterParameters


def order_clusters(sorted_df: pandas.DataFrame, options: OrderClusterParameters) -> Dict[str, List[str]]:
	"""
		Orders genotypes by which background they belong to.
	Parameters
	----------
	sorted_df: pandas.DataFrame
		A dataframe of sorted genotypes based on when the genotype was first detected and first fixed.
	options: OrderClusterParameters
		Parameters for the program to use.

	Returns
	-------
	ClusterType
	"""
	# By default the backgrounds should occupy the first n lines of the dataframe
	initial_background = sorted_df.iloc[0]
	genotype_nests = Cluster(initial_background, timepoints = sorted_df)
	for unnested_label, unnested_trajectory in sorted_df[1:].iterrows():
		logger.info(f"Nesting {unnested_label}")
		genotype_deltas = list()
		# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
		test_table = sorted_df[:unnested_label].iloc[::-1]
		table_of_checks = checks.apply_genotype_checks_to_table(unnested_trajectory, test_table, options)
		table_of_checks = table_of_checks.drop(unnested_label)

		for nested_label, row in table_of_checks.iterrows():
			if nested_label == unnested_label: continue
			additive_check, subtractive_check, delta, area_ratio, area_difference = row
			# logging.info(f"{unnested_label}\t{nested_label}\t{additive_check}\t{subtractive_check}\t{delta}")
			area_check = area_ratio > 0.97
			derivative_check = delta > options.derivative_check_cutoff
			full_check = area_check and derivative_check and additive_check
			logging.info(f"{unnested_label}|{nested_label} Full Check: {full_check}")
			logging.info(f"{unnested_label}|{nested_label} Area Check: {area_check} ({area_ratio})")
			logging.info(f"{unnested_label}|{nested_label} Additive Check: {additive_check}")
			logging.info(f"{unnested_label}|{nested_label} Derivative check: {derivative_check}")
			logging.info(f"{unnested_label}|{nested_label} Area Difference: {area_difference}")
			if area_check and delta > options.derivative_check_cutoff and additive_check:
				# Most likely the background of the current genotype.
				logging.info(f"Complete Check: True")
				genotype_nests.add_genotype_to_background(unnested_label, nested_label, True)
				break
			if area_difference < -.1:
				# The unnested trajectory is larger than the nested trajectory it is being compared against.
				continue
			genotype_deltas.append((nested_label, delta))
			# if subtractive_check:
			#	continue
			if delta > options.derivative_check_cutoff:
				# They are probably on the same background.
				# Need to do one last check: these two genotypes cannot sum to larger than the background.
				genotype_nests.add_genotype_to_background(unnested_label, nested_label)
			elif delta < -options.derivative_check_cutoff:
				# They are anti-correlated.
				continue
			# break
			if area_check:
				# A candidate background
				genotype_nests.add_genotype_to_background(unnested_label, nested_label)
				continue
			if additive_check:# and False:
				# Possible background
				genotype_nests.add_genotype_to_background(unnested_label, nested_label)


		is_member = genotype_nests.is_a_member(unnested_label)
		if not is_member:
			logger.info(f"Not a member: {unnested_label}")
			# if it hasn't been matched with a background, two things can happen
			# (1) IF it is logically possible (i.e., if the sum of it
			# and all existing backgrounds at each time point is ~1 or less),
			# then we can make it its own background. otherwise put
			# it in with the genotype it's most correlated with
			# Find the test that correlated the most with this genotype
			result = background_heuristic(genotype_nests, genotype_deltas, unnested_trajectory, options)
			if result:
				genotype_nests.add_genotype_to_background(unnested_label, result)
	logger.info("The final backgrounds:")
	for k, v in genotype_nests.nests.items():
		vv = "|".join(v)
		logger.info(f"{k}\t{vv}")
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
	back = [unnested_label]
	correlated_label, correlated_delta = get_maximum_genotype_delta(genotype_deltas)
	current_background_total = genotype_nests.get_sum_of_backgrounds()

	try:
		correlated_background = genotype_nests.get(correlated_label)
	except KeyError:
		correlated_background = None

	exceeds_background_limit = (current_background_total > options.new_background_base_cutoff).sum() > 1
	under_cuttoff_limit = (current_background_total < options.new_background_significant_cutoff).sum() == 0

	if exceeds_background_limit or not under_cuttoff_limit:
		# add as background
		result = back
	elif correlated_delta > 0:
		back = correlated_background + back
		result = back
	else:
		# message = 'SOMETHING HAS GONE HORRIBLY WRONG FOR CLUSTER ' + unnested_label
		result = None

	return result[-1]


if __name__ == "__main__":
	pass
