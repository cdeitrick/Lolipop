import pytest
import pandas
from dataio.trajectories import import_table_from_string
from inheritance import timepoint_detection
@pytest.fixture
def transposed_genotypes()->pandas.DataFrame:
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
	table = import_table_from_string(genotype_table_string, index = 'Genotype')
	return table.transpose()

@pytest.fixture
def transposed_mouse_genotypes()->pandas.DataFrame:
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
	t = import_table_from_string(table, index = 'Genotype')
	return t.transpose()

def test_get_first_detected_timepoint(transposed_genotypes, transposed_mouse_genotypes):
	expected = pandas.Series({
		'genotype-1': 25,
		'genotype-2':17,
		'genotype-3': 75,
		'genotype-4': 44,
		'genotype-5': 44,
		'genotype-6': 44,
		'genotype-7': 44,
		'genotype-8': 44,
		'genotype-9': 25,
		'genotype-10': 44,
		'genotype-11': 66,
		'genotype-12': 66,
		'genotype-13':25,
		'genotype-14': 75,
		'genotype-15': 44
	}, name = 'firstDetected')
	expected = expected.astype(str)
	result = timepoint_detection.get_first_detected_timepoint(transposed_genotypes, 0.03)
	assert expected.to_dict() == result.to_dict()

	expected = pandas.Series({
		'genotype-1': 2,
		'genotype-2':1,
		'genotype-3': 1,
		'genotype-4': 0,
		'genotype-5': 5,
		'genotype-6': 4,
		'genotype-7': 6,
		'genotype-8': 7,
		'genotype-9': 1,
		'genotype-10': 4,
		'genotype-11': 1,
		'genotype-12': 0,
		'genotype-13':7,
		'genotype-14': 8,
		'genotype-15': 8
	}, name = 'firstDetected')
	expected = expected.astype(str)
	result = timepoint_detection.get_first_detected_timepoint(transposed_mouse_genotypes, 0.03)
	assert expected.to_dict() == result.to_dict()

def test_get_first_significant_timepoint(transposed_genotypes, transposed_mouse_genotypes):
	expected_mouse = {
		'genotype-1':  3,
		'genotype-2':  1,
		'genotype-3':  3,
		'genotype-4':  0,
		'genotype-5':  5,
		'genotype-6':  5,
		'genotype-7':  9,
		'genotype-8':  9,
		'genotype-9':  10,
		'genotype-10': 10,
		'genotype-11': 10,
		'genotype-12': 1,
		'genotype-13': 10,
		'genotype-14': 10,
		'genotype-15': 10
	}
	expected_mouse = {k:str(v) for k,v in expected_mouse.items()}
	mouse_result = timepoint_detection.get_first_significant_timepoint(transposed_mouse_genotypes, 0.15)

	assert expected_mouse == mouse_result.to_dict()

