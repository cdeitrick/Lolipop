from typing import List, Tuple, Dict

import pandas

from loguru import logger
try:
	from muller.inheritance.checks import apply_genotype_checks
	from muller.inheritance import checks
	from muller.inheritance import scoring
	from muller.inheritance.cluster import Cluster
	from muller.options import OrderClusterParameters
except ModuleNotFoundError:
	from . import scoring
	from . import checks
	from .checks import apply_genotype_checks
	from .cluster import Cluster
	from options import OrderClusterParameters


def order_clusters(sorted_df: pandas.DataFrame, dlimit:float, additive_cutoff:float,derivative_cutoff:float, known_ancestry:Dict[str,str] = None) -> Cluster:
	"""
		Orders genotypes by which background they belong to.
	Parameters
	----------
	sorted_df: pandas.DataFrame
		A dataframe of sorted genotypes based on when the genotype was first detected and first fixed.
	dlimit: float
		The detection limit
	additive_cutoff: float
		Used when testing whther a nested genotype is consistently greater than an unnested genotype
	derivative_cutoff: float
		Used when testing whether two genotypes are correlated, not correlated, or anticorrelated. correlated/anticorrelated genotypes
		must have a covariance outside the range [-`derivative_cutoff`, `derivative_cutoff`].
	known_ancestry: Dict[str,str]
		Manually-assigned ancestry values. For now, the parent genotype is automatically assigned to the root genotype to prevent
		circular links from forming.

	Returns
	-------
	ClusterType
	"""
	# By default the backgrounds should occupy the first n lines of the dataframe
	if known_ancestry is None: known_ancestry = dict()

	initial_background = sorted_df.iloc[0]
	genotype_nests = Cluster(initial_background, timepoints = sorted_df)
	if known_ancestry:
		logger.info(f"Found user-given ancestries.")
	for identity, parent in known_ancestry.items():
		# TODO need to add a way to prevent circular ancestry links when a user manually assigns ancestry. Current workaround forces the manual parent to be in the root background.
		# TODO Also need to add a way to specify a genotype by one of the linked annotations.
		logger.info(f"Adding {parent} as a potential background for {identity}")
		genotype_nests.add_genotype_to_background(parent, 'genotype-0', priority = 100)
		genotype_nests.add_genotype_to_background(identity, parent, priority = 100) # Dummy priority so that it is selected before other backgrounds.
	for unnested_label, unnested_trajectory in sorted_df[1:].iterrows():
		logger.debug(f"Nesting {unnested_label}")
		# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
		test_table = sorted_df[:unnested_label].iloc[::-1]
		for nested_label, nested_genotype in test_table.iterrows():
			if nested_label == unnested_label: continue
			score_additive = scoring.calculate_additive_score(nested_genotype, unnested_trajectory, additive_cutoff)
			score_derivative = scoring.calculate_derivative_score(nested_genotype, unnested_trajectory, detection_cutoff = dlimit, cutoff = derivative_cutoff)
			score_area = scoring.calculate_area_score(nested_genotype, unnested_trajectory)

			total_score = score_additive + score_derivative + score_area
			logger.debug(f"{unnested_label}\t{nested_label}\t{total_score}")
			genotype_nests.add_genotype_to_background(unnested_label, nested_label, total_score)
	logger.debug("Final Ancestry:")
	for genotype_label in sorted_df.index:
		candidate = genotype_nests.get_highest_priority(genotype_label)
		logger.debug(f"{genotype_label}\t{candidate}")
	return genotype_nests


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
	import dataio
	input_string = """
	Genotype	0	1	2	3	4	5	6	7	8	9	10
	genotype-9	0.01	0.279	0.341	0.568	0.708	0.913	0.756	0.455	0.399	0.13	0.041
	genotype-6	0	0	0.045	0.197	0.261	0.096	0.26	0.596	0.66	0.877	0.969
	genotype-18	0	0	0	0	0	0.247	0.388	0.215	0.403	0.141	0.028
	genotype-21	0	0	0	0	0.148	0.384	0.344	0.289	0.333	0.146	0.031
	genotype-17	0	0	0	0	0	0	0.084	0.12	0.124	0.343	0.398
	genotype-20	0	0	0	0	0	0	0	0.077	0.018	0.239	0.308
	genotype-5	0	0	0	0	0	0	0	0	0.236	0.121	0.034
	genotype-16	0.085	0.112	0.16	0	0	0	0	0	0	0	0
	genotype-1	0	0.056	0.101	0.174	0	0	0	0	0	0	0
	genotype-4	0	0	0	0	0	0	0	0.043	0.099	0.146	0.1275
	genotype-11	0.278	0.277	0.224	0.195	0	0	0	0	0	0	0
	genotype-8	0.213	0.265	0.184	0	0	0	0	0	0	0	0
	genotype-12	0.027	0.059	0.0325	0.008	0	0	0	0	0	0	0
	genotype-10	0.013	0	0.026	0.008	0	0	0	0	0	0	0.041
	genotype-7	0	0.088	0.036	0.046	0	0.059	0.052	0	0.073	0	0
	genotype-14	0	0.021	0.026	0.02	0.031	0	0	0.032	0.035	0.016	0
	genotype-13	0	0	0	0	0.072	0.047	0.057	0	0	0	0
	genotype-19	0	0	0	0	0	0	0.03	0.085	0.037	0.019	0
	genotype-3	0	0	0	0	0	0	0.007666666666667	0.006333333333333	0.033	0.043333333333333	0.026333333333333
	genotype-15	0	0	0	0	0	0	0.007	0.038	0.023666666666667	0.023333333333333	0.041
	genotype-2	0	0	0	0	0	0	0	0	0.067	0.068	0.074
	"""
	table = dataio.import_table(input_string, index = 'Genotype')
	order_clusters(table, additive_cutoff = 0.03, derivative_cutoff = 0.03, dlimit = 0.03)
