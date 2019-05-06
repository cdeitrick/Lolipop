import pandas
import pytest

from muller.clustering import filters
from muller.dataio import import_table


@pytest.fixture
def psdata() -> pandas.DataFrame:
	string = """
		Trajectory	3	4	6	7	9	10	12
		1	0	0	0	0.286	0.281	0.054	0.167	0
		2	0	0	0	0.062	0.215	0	0.074	0
		3	0.198	0.758	0.696	0.084	0.065	0.204	0	0
		4	0	0	0.069	0.282	0.168	0	0	0
		5	0	0	0.095	0	0.117	0.416	0.36	0
		6	0	0.078	0.118	0	0	0	0	0
		7	0	0	0	0	0	0.093	0	0
		8	0	0	0	0	0.227	0	0	0
		9	0	0	0	0	0.287	0	0	0
		10	0	0	0	0	0	0	0.132	0
		11	0	0	0	0.171	0	0	0	0
		12	0	0	1	0	0	0	0	0
		13	0	0	0	0	0	0	0.238	0
		14	0	0	0	0	0	0	0.055	0
		15	0.059	0	0.053	0.159	0.18	0.067	0.12	0
		16	0.059	0	0	0.103	0.087	0.465	0.294	0
		17	0.194	0	0	0.207	0.362	0.135	0.216	0
		18	0	0	0	0.128	0.114	0	0.125	0
		19	0	0	0	0.22	0.29	0.055	0.113	0
		20	0	0.06	0	0	0	0	0	0
	"""
	df = import_table(string, index = "Trajectory")
	return df


@pytest.fixture
def genotype_table_a() -> pandas.DataFrame:
	string = """Genotype	0	3	4	6	7	9	10	12
		genotype-1	0	0.198	0.758	0.696	0.084	0.065	0.204	0
		genotype-2	0	0	0	1	0	0	0	0
		genotype-3	0	0.0295	0	0.0475	0.0515	0.102	0.4405	0.327
		genotype-4	0	0	0	0	0	0	0	0.141666666666667
		genotype-5	0	0	0.046	0.039333333333333	0	0	0.031	0
		genotype-6	0	0.06325	0	0.01325	0.218	0.27825	0.07775	0.154
		genotype-7	0	0	0	0.0115	0.107166666666667	0.1685	0	0.033166666666667
	"""
	data = import_table(string, index = "Genotype")
	data.columns = [int(i) for i in data.columns]
	return data


@pytest.fixture
def genotype_table_b() -> pandas.DataFrame:
	string = """                                                                 
		Genotype	0	3	4	6	7	9	10	12
		genotype-1	0.0	0.19800	0.758	0.696000	0.084000	0.06500	0.20400	0.000000
		genotype-2	0.0	0.02950	0.000	0.047500	0.051500	0.10200	0.44050	0.327000
		genotype-3	0.0	0.06325	0.000	0.013250	0.218000	0.27825	0.07775	0.154000
		genotype-4	0.0	0.00000	0.023	0.019667	0.000000	0.00000	0.01550	0.070833
		genotype-5	0.0	0.00000	0.000	0.011500	0.107167	0.16850	0.00000	0.033167
	"""
	data = import_table(string, index = 'Genotype')
	data.columns = [int(i) for i in data.columns]
	return data


def test_trajectory_filters(psdata):
	dlimit = 0.03
	flimit = 0.97

	result = filters.filter_trajectories(psdata, dlimit, flimit)
	expected_index = list(map(str, range(1, 21)))
	expected_index.remove("12")
	assert expected_index == list(result.index)


def test_get_fuzzy_backgrounds(genotype_table_a):
	test_background, flimit = filters.get_fuzzy_backgrounds(genotype_table_a, [1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0])

	assert pytest.approx(flimit, 0.7)
	assert len(test_background) == 1
	assert list(test_background.index) == ['genotype-1']


def test_find_first_invalid_genotype_a(genotype_table_a):
	backgrounds = genotype_table_a.loc['genotype-1'].to_frame().transpose()
	dlimit, flimit = 0.03, 0.7

	first_invalid_genotype = filters.find_first_invalid_genotype(genotype_table_a, backgrounds, dlimit, flimit, False)

	assert first_invalid_genotype is None


def test_find_first_invalid_genotype_b(genotype_table_b):
	backgrounds = genotype_table_b.loc['genotype-1'].to_frame().transpose()
	first_invalid_genotype = filters.find_first_invalid_genotype(genotype_table_b, backgrounds, 0.03, 0.7, False)

	assert first_invalid_genotype == None
