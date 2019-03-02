import pytest
from dataio import import_table
from clustering.average import _calculate_mean_frequencies_of_trajectories
import pandas

@pytest.fixture
def trajectory_table():
	trajectory_table_string = """
		Trajectory	0	17	25	44	66	75	90	genotype
		1	0	0	0.261	1	1	1	1	genotype-1
		20	0	0	0	0.138	0.295	0	0.081	genotype-10
		4	0	0	0	0	0.211	0.811	0.813	genotype-11
		8	0	0	0	0	0.345	0.833	0.793	genotype-11
		16	0	0	0	0	0.209	0.209	0	genotype-12
		13	0	0	0	0	0.258	0.057	0.075	genotype-12
		15	0	0	0.066	0.104	0.062	0	0	genotype-13
		11	0	0	0	0.108	0.151	0	0	genotype-13
		17	0	0	0	0	0	0.266	0.312	genotype-14
		9	0	0	0	0	0	0.269	0.34	genotype-14
		18	0	0	0	0.115	0	0.131	0	genotype-15
		21	0	0	0	0.114	0	0.11	0.123	genotype-15
		14	0	0.38	0.432	0	0	0	0	genotype-2
		6	0	0	0	0	0	1	1	genotype-3
		2	0	0	0	0.525	0.454	0.911	0.91	genotype-4
		3	0	0	0	0.147	0.45	0.924	0.887	genotype-5
		7	0	0	0	0.273	0.781	1	1	genotype-6
		19	0	0	0	0.188	0.171	0.232	0.244	genotype-7
		5	0	0	0	0.403	0.489	0.057	0.08	genotype-8
		10	0	0	0.117	0	0	0	0.103	genotype-9
		12	0	0.125	0	0.153	0.181	0.175	0.191	filtered
	"""
	return import_table(trajectory_table_string, index = 'Trajectory')

@pytest.fixture
def genotype_table():
	genotype_table_string = """
		Genotype	0	17	25	44	66	75	90
		genotype-1	0	0	0.261	1	1	1	1
		genotype-2	0	0.38	0.432	0	0	0	0
		genotype-3	0	0	0	0	0	1	1
		genotype-4	0	0	0	0.525	0.454	0.911	0.91
		genotype-5	0	0	0	0.147	0.45	0.924	0.887
		genotype-6	0	0	0	0.273	0.781	1	1
		genotype-7	0	0	0	0.188	0.171	0.232	0.244
		genotype-8	0	0	0	0.403	0.489	0.057	0.08
		genotype-9	0	0	0.117	0	0	0	0.103
		genotype-10	0	0	0	0.138	0.295	0	0.081
		genotype-11	0	0	0	0	0.278	0.822	0.803
		genotype-12	0	0	0	0	0.2335	0.133	0.0375
		genotype-13	0	0	0.033	0.106	0.1065	0	0
		genotype-14	0	0	0	0	0	0.2675	0.326
		genotype-15	0	0	0	0.1145	0	0.1205	0.0615
	"""
	table = import_table(genotype_table_string, index = 'Genotype')
	table = table.astype(float)
	return table
def test_calculate_mean_genotype(genotype_table, trajectory_table):
	groups = trajectory_table.groupby(by = 'genotype')
	for genotype_label, genotype_group in groups:
		if genotype_label == 'filtered': continue
		truth_genotype = genotype_table.loc[genotype_label]
		mean_genotype = _calculate_mean_frequencies_of_trajectories(genotype_label, genotype_group, ['a', 'b', 'c'])
		members = mean_genotype.pop('members')
		mean_genotype = mean_genotype.astype(float)
		assert 'a|b|c' == members
		pandas.testing.assert_series_equal(truth_genotype, mean_genotype)