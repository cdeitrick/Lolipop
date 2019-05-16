import pandas
import pytest

from muller.clustering import filters
from muller.dataio import import_table


@pytest.fixture
def genotypes() -> pandas.DataFrame:
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
	t = import_table(genotype_table_string, index = 'Genotype')
	t = t.astype(float)
	return t


@pytest.fixture
def trajectory_filter() -> filters.TrajectoryFilter:
	f = filters.TrajectoryFilter(detection_cutoff = 0.03, fixed_cutoff = 0.97)
	return f


@pytest.fixture
def genotype_filter() -> filters.GenotypeFilter:
	g = filters.GenotypeFilter(detection_cutoff = 0.03, fixed_cutoff = 0.97, frequencies = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0])
	return g


def test_remove_single_point_series(trajectory_filter):
	table = pandas.DataFrame([
		[0, 0, 0, 1, 0, 0],
		[1, 3, 2, 0, 1, 1],
		[0.03, 1, 0, 0, 0, 0]
	])
	single_point_series = trajectory_filter._remove_single_point_background(table)

	assert [0] == list(single_point_series)

	trajectory_filter.dlimit = 0.04
	single_point_series = trajectory_filter._remove_single_point_background(table)

	assert [0, 2] == list(single_point_series)


def test_get_first_timepoint_above_cutoff(genotype_filter):
	series = pandas.Series([0, .03, .04, .1, .2, .3, .4, .5, .6, .7, .8])

	assert 2 == genotype_filter.get_first_timepoint_above_cutoff(series, 0.03)
	assert 4 == genotype_filter.get_first_timepoint_above_cutoff(series, 0.1)
	assert 1 == genotype_filter.get_first_timepoint_above_cutoff(series, 0)


def test_get_fuzzy_backgrounds(genotypes, genotype_filter):
	expected = """
		Genotype	0	17	25	44	66	75	90
		genotype-1	0	0	0.261	1	1	1	1
		genotype-4	0	0	0	0.525	0.454	0.911	0.91
		genotype-5	0	0	0	0.147	0.45	0.924	0.887
		genotype-6	0	0	0	0.273	0.781	1	1"""
	expected_table = import_table(expected, index = 'Genotype')
	genotype_filter.frequencies = [0.9]
	backgrounds = genotype_filter.get_fuzzy_backgrounds(genotypes)
	expected_table = expected_table.astype(float)
	assert pytest.approx(genotype_filter.fuzzy_fixed_cutoff == 0.9)

	pandas.testing.assert_frame_equal(expected_table, backgrounds)
