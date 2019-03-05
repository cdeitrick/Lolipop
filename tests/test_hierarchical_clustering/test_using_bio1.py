import pytest
import pandas
from dataio import import_table
from clustering.methods.hierarchical_method import hierarchical_method
from clustering.metrics import calculate_pairwise_metric, PairwiseCalculationCache
@pytest.fixture
def bio1() -> pandas.DataFrame:
	string = """
		Trajectory	0	3	4	6	7	9	10	12
		1	0	0.653	1	1	1	0.91	0.907	1
		2	0	0	0.646	0.777	0.89	0.512	0.135	0.546
		3	0	0.054	0.19	0.131	0	0	0.053	0.124
		4	0	0.08	0	0	0	0	0	0.226
		5	0	0.187	0	0	0	0	0.117	0
		6	0	0	0	0	0	0	0.759	0
	"""
	return import_table(string, index = 'Trajectory')

def test_clustering(bio1):
	expected = [
		["1"],
		["2"],
		["3"],
		["4", "5"],
		["6"]
	]
	pair_array = calculate_pairwise_metric(bio1, .03, 0.97, 'binomial')
	_, result = hierarchical_method(PairwiseCalculationCache(pair_array), similarity_cutoff = 0.03)

	assert expected == result

