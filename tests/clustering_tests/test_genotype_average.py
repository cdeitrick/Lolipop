from pathlib import Path

import pandas
import pytest
from loguru import logger
from muller.clustering.generate_genotypes import ClusterMutations
from muller.dataio import import_table
from tests import filenames

DATA_FOLDER = Path(__file__).parent.parent / "data" / "tables"


@pytest.fixture
def genotype_generator() -> ClusterMutations:
	generator = ClusterMutations(
		metric = 'binomial',
		dlimit = 0.03,
		flimit = 0.97,
		starting_genotypes = [],
	)
	return generator


@pytest.mark.parametrize("filename", list(filenames.generic_tables.values()) + list(filenames.model_tables.values()))
def test_calculate_mean_genotype(genotype_generator, filename):
	trajectories = import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	genotypes = import_table(filename, sheet_name = 'genotype', index = 'Genotype')

	trajectories['Genotype'] = ['genotype-' + i.split('-')[1] for i in trajectories.index]
	logger.debug(trajectories['Genotype'])
	groups = trajectories.groupby(by = "Genotype")
	for genotype_label, group in groups:
		expected = genotypes.loc[genotype_label].astype(float)
		mean_genotype = genotype_generator._calculate_mean_frequencies_of_trajectories(genotype_label, group,
			group.index)
		del mean_genotype['members']
		mean_genotype = mean_genotype.astype(float)

		pandas.testing.assert_series_equal(mean_genotype, expected, check_index_type = False)
