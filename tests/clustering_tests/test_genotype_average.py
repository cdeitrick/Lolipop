from pathlib import Path

import pandas
import pytest

from muller.clustering.generate_genotypes import ClusterMutations
from muller.dataio import import_table

DATA_FOLDER = Path(__file__).parent.parent / "data" / "tables"

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


@pytest.fixture
def genotype_generator() -> ClusterMutations:
	generator = ClusterMutations(
		metric = 'binomial',
		dlimit = 0.03,
		slimit = 0.15,
		flimit = 0.97,
		pvalue = 0.05,
		starting_genotypes = [],
	)
	return generator


@pytest.mark.parametrize("filename", generic_tables + model_tables)
def test_calculate_mean_genotype(genotype_generator, filename):
	trajectories = import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	genotypes = import_table(filename, sheet_name = 'genotype', index = 'Genotype')

	trajectories['Genotype'] = ['genotype-' + i.split('-')[1] for i in trajectories.index]

	groups = trajectories.groupby(by = "Genotype")
	for genotype_label, group in groups:
		expected = genotypes.loc[genotype_label].astype(float)
		mean_genotype = genotype_generator._calculate_mean_frequencies_of_trajectories(genotype_label, group, group.index)
		del mean_genotype['members']
		mean_genotype = mean_genotype.astype(float)

		pandas.testing.assert_series_equal(mean_genotype, expected, check_index_type = False)
