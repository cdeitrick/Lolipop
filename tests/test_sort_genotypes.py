import pandas
import pytest

from muller.dataio import import_table
from muller.inheritance.genotype_reorder import SortGenotypeTableWorkflow


@pytest.fixture
def genotype_sorter() -> SortGenotypeTableWorkflow:
	return SortGenotypeTableWorkflow(0.03, 0.15, 0.97, [1, 0.97, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0])


@pytest.fixture
def table() -> pandas.DataFrame:
	trajectory_table = """
		Trajectory	0	17	25	44	66	75	90
		1	0	0	0.261	1	1	1	1
		20	0	0	0	0.138	0.295	0	0.081
		16	0	0	0	0	0.209	0.209	0
		13	0	0	0	0	0.258	0.057	0.075
		15	0	0	0.066	0.104	0.062	0	0
		11	0	0	0	0.108	0.151	0	0
		17	0	0	0	0	0	0.266	0.312
		9	0	0	0	0	0	0.269	0.34
		14	0	0.38	0.432	0	0	0	0
		6	0	0	0	0	0	1	1
		2	0	0	0	0.525	0.454	0.911	0.91
		3	0	0	0	0.147	0.45	0.924	0.887
		7	0	0	0	0.273	0.781	1	1
		10	0	0	0.117	0	0	0	0.103
	"""
	t = import_table(trajectory_table, index = 'Trajectory')
	return t


@pytest.fixture
def mouse_table() -> pandas.DataFrame:
	table = """
	Genotype	0	1	2	3	4	5	6	7	8	9	10
	genotype-1	0	0	0.045	0.197	0.261	0.096	0.26	0.596	0.66	0.877	0.969
	genotype-2	0.01	0.279	0.341	0.568	0.708	0.913	0.756	0.455	0.399	0.13	0.041
	genotype-3	0	0.056	0.101	0.174	0	0	0	0	0	0	0
	genotype-4	0.278	0.277	0.224	0.195	0	0	0	0	0	0	0
	genotype-5	0	0	0	0	0	0.247	0.388	0.215	0.403	0.141	0.028
	genotype-6	0	0	0	0	0.148	0.384	0.344	0.289	0.333	0.146	0.031
	genotype-7	0	0	0	0	0	0	0.084	0.12	0.124	0.343	0.398
	genotype-8	0	0	0	0	0	0	0	0.077	0.018	0.239	0.308
	genotype-9	0	0.088	0.036	0.046	0	0.059	0.052	0	0.073	0	0
	genotype-10	0	0	0	0	0.072	0.047	0.057	0	0	0	0
	genotype-11	0.027	0.059	0.0325	0.008	0	0	0	0	0	0	0
	genotype-12	0.149	0.1885	0.172	0	0	0	0	0	0	0	0
	genotype-13	0	0.00525	0.0065	0.005	0.00775	0	0.01275	0.051	0.032	0.0195	0.02175
	genotype-14	0	0	0	0	0	0	0	0.0172	0.1156	0.112	0.0948
	genotype-15	0.001857	0	0.003714	0.001143	0	0	0.003286	0.006571	0.034	0.040286	0.038143
	"""
	t = import_table(table, index = 'Genotype')
	return t


@pytest.fixture
def smalltable() -> pandas.DataFrame:
	trajectory_table = """
		Trajectory	0	17	25	44	66	75	90
		1	0	0	0.261	1	1	1	1
		20	0	0	0	0.138	0.295	0	0.081
		16	0	0	0	0	0.209	0.209	0
		15	0	0	0.066	0.104	0.062	0	0
		9	0	0	0	0	0	0.269	0.34
		10	0	0	0.117	0	0	0	0.103
	"""
	t = import_table(trajectory_table, index = 'Trajectory')
	return t


@pytest.fixture
def genotype_table() -> pandas.DataFrame:
	trajectory_table = """
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
	t = import_table(trajectory_table, index = 'Genotype')
	return t


def test_sort_genotypes(table, genotype_sorter):
	expected = """
		Trajectory	0	17	25	44	66	75	90
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
	expected_result = import_table(expected, index = 'Trajectory')

	result = genotype_sorter.run(table)
	expected_result.index.name = None
	pandas.testing.assert_frame_equal(result, expected_result)


def test_sort_genotypes_class(table, genotype_sorter):
	expected = """
		Trajectory	0	17	25	44	66	75	90
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
	expected_result = import_table(expected, index = 'Trajectory')
	result = genotype_sorter.run(table)
	expected_result.index.name = None
	pandas.testing.assert_frame_equal(result, expected_result)


def test_sort_genotypes_with_initial_values(mouse_table, genotype_sorter):
	expected = """
		Genotype	0	1	2	3	4	5	6	7	8	9	10
		genotype-1	0	0	0.045	0.197	0.261	0.096	0.26	0.596	0.66	0.877	0.969
		genotype-2	0.01	0.279	0.341	0.568	0.708	0.913	0.756	0.455	0.399	0.13	0.041
		genotype-5	0	0	0	0	0	0.247	0.388	0.215	0.403	0.141	0.028
		genotype-8	0	0	0	0	0	0	0	0.077	0.018	0.239	0.308
		genotype-6	0	0	0	0	0.148	0.384	0.344	0.289	0.333	0.146	0.031
		genotype-7	0	0	0	0	0	0	0.084	0.12	0.124	0.343	0.398
		genotype-4	0.278	0.277	0.224	0.195	0	0	0	0	0	0	0
		genotype-12	0.149	0.1885	0.172	0	0	0	0	0	0	0	0
		genotype-3	0	0.056	0.101	0.174	0	0	0	0	0	0	0
		genotype-14	0	0	0	0	0	0	0	0.0172	0.1156	0.112	0.0948
		genotype-11	0.027	0.059	0.0325	0.008	0	0	0	0	0	0	0
		genotype-15	0.001857	0	0.003714	0.001143	0	0	0.003286	0.006571	0.034	0.040286	0.038143
		genotype-9	0	0.088	0.036	0.046	0	0.059	0.052	0	0.073	0	0
		genotype-13	0	0.00525	0.0065	0.005	0.00775	0	0.01275	0.051	0.032	0.0195	0.02175
		genotype-10	0	0	0	0	0.072	0.047	0.057	0	0	0	0
	"""
	expected_result = import_table(expected, index = 'Genotype')
	expected_result.index.name = None
	from loguru import logger
	logger.info(f"{genotype_sorter.breakpoints}")
	result = genotype_sorter.run(mouse_table)
	pandas.testing.assert_frame_equal(expected_result, result)
