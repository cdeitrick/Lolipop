from pathlib import Path
import pytest
from muller import dataio
import pandas
from loguru import logger
from muller.clustering import ClusterMutations

@pytest.fixture
def cluster()->ClusterMutations:
	c = ClusterMutations(
		method = 'hierarchy',
		metric = 'binomial',
		dlimit = 0.03,
		flimit = 0.97,
		pvalue = 0.05,
		dbreakpoint = 0.15, #won't be used since all tests use hierarchical method.
		breakpoints = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
	)

	return c

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

def helper_get_expected_members(table:pandas.DataFrame):
	# Groups each trajectory into the expected genotype.
	table = table.reset_index() # Because we ned to acces the 'Trajectory' column
	if 'Genotype' not in table.columns:
		table['Genotype'] = ['genotype-'+i.split('-')[1] for i in table['Trajectory'].values]
	groups = table.groupby(by = 'Genotype')
	members = dict()
	for genotype_label, group in groups:
		members[genotype_label] = list(group['Trajectory'])
	members['genotype-filtered'] = []
	return members

@pytest.mark.parametrize("filename", generic_tables)
def test_clustering_algorithm_on_generic_tables(cluster, filename):
	cluster.pvalue = 0.05
	trajectories = dataio.import_table(filename, sheet_name = 'trajectory', index = 'Trajectory')
	genotypes = dataio.import_table(filename, sheet_name = 'genotype', index = 'Genotype')
	expected_members = helper_get_expected_members(trajectories)

	if 'Genotype' in trajectories.columns:
		del trajectories['Genotype'] # Make sure the table only as numeric values.

	mean, members = cluster.run(trajectories)

	assert sorted(members.values()) == sorted(expected_members.values())

