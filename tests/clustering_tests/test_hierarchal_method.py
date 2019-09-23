from pathlib import Path

import pandas
import pytest
from loguru import logger

from muller import dataio
from muller.clustering.methods import HierarchalCluster, format_linkage_matrix
from muller.clustering.metrics import DistanceCache, DistanceCalculator
from .. import filenames
from typing import Dict, List

@pytest.fixture
def calculator()->DistanceCalculator:
	return DistanceCalculator(0.03, 0.97, 'binomial')
@pytest.fixture
def hierarchy()->HierarchalCluster:
	return HierarchalCluster()


def helper_get_expected_genotypes(trajectories: pandas.Series)->Dict[str,List[str]]:
	expected = dict()
	for value in trajectories:
		genotype_name = value.split('-')[1]
		expected[genotype_name] = expected.get(genotype_name, []) + [value]

	return expected


@pytest.mark.parametrize("filename", list(filenames.generic_tables_with_trajectories.values()) + list(filenames.real_tables.values()))
def test_default_hierarchal_method(calculator, hierarchy, filename):
	trajectories = dataio.import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	distances = calculator.run(trajectories)
	expected = helper_get_expected_genotypes(trajectories.index)
	pair_array = DistanceCache(distances)

	result, Z = hierarchy.run(pair_array)
	assert sorted(result) == sorted(expected.values())

def test_format_linkage_matrix():
	string = """
			left	right	distance	observations
		7	18	0.033999578315858	2
		13	17	0.17508789172405	2
		8	11	0.199140037464566	2
		2	5	0.238709275657774	2
		10	3	0.278982267870099	2
		9	12	0.370131434040108	2
		23	6	0.528725037646305	3
		22	21	0.624301321943297	4
		26	1	0.708258027601832	5
		24	16	0.760211707999897	3
		14	25	0.785856622308224	4
		15	20	0.987877097254146	3
		29	27	1.09399939089444	9
		31	19	1.35849899883873	11
		30	28	1.36249792376227	6
		4	32	1.49918800504887	12
		33	0	2.12548563025534	7
		34	35	4.94288036257192	19
	"""

	expected_table = """
		left	right	distance	observations	clusterId
		7	18	0.033999578315858	2	19
		13	17	0.17508789172405	2	20
		8	11	0.199140037464566	2	21
		2	5	0.238709275657774	2	22
		10	3	0.278982267870099	2	23
		9	12	0.370131434040108	2	24
		23	6	0.528725037646305	3	25
		22	21	0.624301321943297	4	26
		26	1	0.708258027601832	5	27
		24	16	0.760211707999897	3	28
		14	25	0.785856622308224	4	29
		15	20	0.987877097254146	3	30
		29	27	1.09399939089444	9	31
		31	19	1.35849899883873	11	32
		30	28	1.36249792376227	6	33
		4	32	1.49918800504887	12	34
		33	0	2.12548563025534	7	35
		34	35	4.94288036257192	19	36
		"""
	test_table = dataio.import_table(string)
	expected_table = dataio.import_table(expected_table)
	test_table = format_linkage_matrix(test_table, 19)

	# Make sure the test table's columns are aligned the same way.

	pandas.testing.assert_frame_equal(test_table, expected_table)