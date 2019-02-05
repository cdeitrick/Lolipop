import pandas
import pytest

from import_data import import_table_from_string
from inheritance.sort_genotypes import _get_timepoint_above_threshold, sort_genotypes


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
	t = import_table_from_string(trajectory_table, index = 'Trajectory')
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
	t = import_table_from_string(trajectory_table, index = 'Trajectory')
	return t

@pytest.fixture
def genotype_table()->pandas.DataFrame:
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
	t = import_table_from_string(trajectory_table, index = 'Genotype')
	return t


@pytest.mark.parametrize(
	"threshold,expected",
	[
		(0.03, {'1': '25', '20': '44', '16': '66', '15': '25', '9': '75', '10': '25'}),
		(0.97, {'1': '44', '20': '0', '16': '0', '15': '0', '9': '0', '10': '0'})
	]
)
def test_timepoint_above_threshold(smalltable, threshold, expected):
	result = _get_timepoint_above_threshold(smalltable.T, threshold)
	assert expected == result.to_dict()

def test_sort_genotypes(table):
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
	expected_result = import_table_from_string(expected, index = 'Trajectory')
	class TestOptions:
		def __init__(self):
			self.detection_breakpoint = 0.03
			self.fixed_breakpoint = 0.97
			self.significant_breakpoint = 0.15
			self.frequency_breakpoints = [1,.9,.8,.7,.6,.5,.4,.3,.2,.1,0]
	result = sort_genotypes(table, TestOptions())
	expected_result.index.name = None
	pandas.testing.assert_frame_equal(result, expected_result)
