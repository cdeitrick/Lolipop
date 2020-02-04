from pathlib import Path

import pandas
import pytest

from muller import dataio
from muller.clustering import ClusterMutations
from .. import filenames


@pytest.fixture
def cluster() -> ClusterMutations:
	value = ClusterMutations(
		metric = 'binomial',
		dlimit = 0.03,
		slimit = 0.15,
		flimit = 0.97,
		pvalue = 0.05
	)

	return value


DATA_FOLDER = Path(__file__).parent.parent / "data" / "tables"


def helper_get_expected_members(table: pandas.DataFrame):
	# Groups each trajectory into the expected genotype.
	members = dict()
	for element in table.index:
		genotype_name = element.split('-')[1]
		members[genotype_name] = members.get(genotype_name, []) + [element]
	return members


@pytest.mark.parametrize("filename", filenames.generic_tables_with_trajectories.values())
def test_clustering_algorithm_on_generic_tables(cluster, filename):

	trajectories = dataio.import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	if len(trajectories) == 10:
		# The table with 10 genotypes shouldn't be used as a trajectory table.
		return None

	expected_members = helper_get_expected_members(trajectories)

	result = cluster.run(trajectories, distance_cutoff = 0.2)

	assert sorted(result.genotype_members.values()) == sorted(expected_members.values())


@pytest.mark.parametrize("filename", filenames.real_tables.values())
def test_clustering_algorithm_on_real_tables(cluster, filename):
	trajectories = dataio.import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	expected_members = helper_get_expected_members(trajectories)

	result = cluster.run(trajectories, distance_cutoff = 0.2)

	assert sorted(result.genotype_members.values()) == sorted(expected_members.values())
