import itertools
from typing import Dict, Tuple

import pandas
import pytest

from dataio.trajectories import import_table
from muller.clustering.metrics import distance
from muller.clustering.metrics.pairwise_calculation_cache import PairwiseCalculationCache
from muller.clustering.methods.hierarchical_method import hierarchical_method

@pytest.fixture
def five_genotypes() -> pandas.DataFrame:
	string = """
	Trajectory	0	1	2	3	5	8	13	21	34
	trajectory-A1	0	0	0	0.1	0.5	0.5	0.61	0.7	0.8
	trajectory-A2	0	0	0	0.06	0.35	0.4	0.65	0.75	0.8
	trajectory-A3	0	0	0	0	0.45	0.5	0.55	0.7	0.85
	trajectory-B1	0	0.1	0.15	0.03	0	0	0	0	0
	trajectory-B2	0	0.07	0.1	0.02	0.01	0	0	0	0
	trajectory-C1	0	0	0	0.3	0.7	1	1	1	1
	trajectory-D1	0	0	0	0	0	0	0.1	0	0
	trajectory-D2	0	0	0	0	0	0	0.07	0	0.01
	trajectory-E1	0	0	0	0	0	0	0.1	0.5	0.5
	trajectory-E2	0	0	0	0	0	0	0.05	0.55	0.5
	trajectory-E3	0	0	0	0	0	0	0	0.6	0.5
	"""
	table = import_table(string, index = 'Trajectory')

	return table


def test_using_minkowski_simple_case(five_genotypes):
	combinations = itertools.combinations(five_genotypes.index, 2)
	minkowski_distances: Dict[Tuple[str, str], float] = dict()
	for left_label, right_label in combinations:
		left_trajectory = five_genotypes.loc[left_label]
		right_trajectory = five_genotypes.loc[right_label]

		d = distance.minkowski_distance(left_trajectory, right_trajectory, 3)
		minkowski_distances[left_label, right_label] = d
		minkowski_distances[right_label, left_label] = d
	cache = PairwiseCalculationCache(minkowski_distances)

	clusters, linkage_matrix = hierarchical_method(cache, 0.05)
	print(len(clusters))
	print(clusters)



