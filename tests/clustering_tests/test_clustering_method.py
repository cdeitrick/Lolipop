from pathlib import Path
import pytest
from muller import dataio
import pandas
from loguru import logger
from muller.clustering import ClusterMutations

@pytest.fixture
def cluster()->ClusterMutations:
	value =  ClusterMutations(
		metric = 'binomial',
		dlimit = 0.03,
		slimit = 0.15,
		flimit = 0.97,
		pvalue = 0.05
	)

	return value
from ..filenames import generic_tables_with_trajectories
DATA_FOLDER = Path(__file__).parent.parent / "data" / "tables"

def helper_get_expected_members(table:pandas.DataFrame):
	# Groups each trajectory into the expected genotype.
	members = dict()
	for element in table.index:
		genotype_name = element.split('-')[1]
		members[genotype_name] = members.get(genotype_name, []) + [element]
	#members['filtered'] = []
	return members

@pytest.mark.parametrize("filename", generic_tables_with_trajectories.values())
def test_clustering_algorithm_on_generic_tables(cluster, filename):
	cluster.pvalue = 0.05
	trajectories = dataio.import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	expected_members = helper_get_expected_members(trajectories)

	result = cluster.run(trajectories, distance_cutoff = None)

	assert sorted(result.genotype_members.values()) == sorted(expected_members.values())

