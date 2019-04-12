import pandas
import pytest

from clustering.metrics.calculation import calculate_pairwise_metric
from dataio import import_table


@pytest.fixture
def b1_data() -> pandas.DataFrame:
	string = """
		Trajectory	0	17	25	44	66	75	90
		1	0.000	0.000	0.261	1.000	1.000	1.000	1.000
		2	0.000	0.000	0.000	0.525	0.454	0.911	0.910
		3	0.000	0.000	0.000	0.147	0.450	0.924	0.887
		4	0.000	0.000	0.000	0.000	0.211	0.811	0.813
		5	0.000	0.000	0.000	0.403	0.489	0.057	0.080
		6	0.000	0.000	0.000	0.000	0.000	1.000	1.000
		7	0.000	0.000	0.000	0.273	0.781	1.000	1.000
		8	0.000	0.000	0.000	0.000	0.345	0.833	0.793
		9	0.000	0.000	0.000	0.000	0.000	0.269	0.340
		10	0.000	0.000	0.117	0.000	0.000	0.000	0.103
		11	0.000	0.000	0.000	0.108	0.151	0.000	0.000
		12	0.000	0.125	0.000	0.153	0.181	0.175	0.191
		13	0.000	0.000	0.000	0.000	0.258	0.057	0.075
		14	0.000	0.380	0.432	0.000	0.000	0.000	0.000
		15	0.000	0.000	0.066	0.104	0.062	0.000	0.000
		16	0.000	0.000	0.000	0.000	0.209	0.209	0.000
		17	0.000	0.000	0.000	0.000	0.000	0.266	0.312
		18	0.000	0.000	0.000	0.115	0.000	0.131	0.000
		19	0.000	0.000	0.000	0.188	0.171	0.232	0.244
		20	0.000	0.000	0.000	0.138	0.295	0.000	0.081
		21	0.000	0.000	0.000	0.114	0.000	0.110	0.123
	"""
	return import_table(string, index = 'Trajectory')


@pytest.fixture
def pairwise_values(b1_data):
	pairwise_values = calculate_pairwise_metric(b1_data, 0.03, 0.97, 'binomial')
	return pairwise_values


@pytest.mark.parametrize("name,closest",
	[
		("4", "8"),
		("2", "3"),
		("17", "9"),
		("5", "20")
	]
)
def test_binomial_distance(pairwise_values, name: str, closest: str):
	candidates = (i for i in pairwise_values.items() if name in i[0])
	pair, value = min(candidates, key = lambda s: s[1])
	assert name in pair and closest in pair
