from pathlib import Path

import pandas
import pytest
from loguru import logger

from muller import dataio
from muller.clustering.methods import hierarchical_method
from muller.clustering.metrics import DistanceCache

DATA_FOLDER = Path(__file__).parent / "data" / "tables"

generic_tables = [
	DATA_FOLDER / "generic.coexistinglineages.xlsx",
	DATA_FOLDER / "generic.genotypes.3.xlsx",
	DATA_FOLDER / "generic.genotypes.5.xlsx",
	DATA_FOLDER / "generic.genotypes.10.xlsx",
	DATA_FOLDER / "generic.small.xlsx"
]

model_tables = [
	DATA_FOLDER / "model.clonalinterferance.xlsx",
	DATA_FOLDER / "model.periodicselection.xlsx",
	DATA_FOLDER / "model.strongselection.xlsx"
]

real_tables = [
	DATA_FOLDER / "real.nature12344-s2.BYB1-G07.xlsx"
]


def helper_run_hierarchal_method(pair_array: DistanceCache):
	similarity_cutoff = 0.05
	cluster_method = 'distance'

	return hierarchical_method(pair_array, similarity_cutoff, cluster_method)


def helper_get_expected_genotypes(trajectories: pandas.DataFrame):
	trajectories['Genotype'] = [i[:-1] for i in trajectories['Trajectory']]
	groups = trajectories.groupby(by = 'Genotype')
	expected = dict()
	for i, g in groups:
		expected[i] = sorted(g['Trajectory'])
	return expected


@pytest.mark.parametrize("filename", [generic_tables[1]])
def test_default_hierarchal_method(filename):
	square_table = dataio.import_table(filename, sheet_name = 'distance', index = 'name')

	trajectories = dataio.import_table(filename, sheet_name = 'trajectory', index = 'Trajectory').reset_index()

	expected = helper_get_expected_genotypes(trajectories)

	logger.info(square_table.columns)
	pair_array = DistanceCache.from_squareform(square_table)

	result, Z = helper_run_hierarchal_method(pair_array)


	assert sorted(result) == sorted(expected.values())
