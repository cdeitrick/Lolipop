import pytest
from pathlib import Path
from muller import dataio, inheritance
DATA_FOLDER = Path(__file__).parent / "data" / "tables_input_genotypes"
from loguru import logger
@pytest.fixture
def lineage_workflow()->inheritance.LineageWorkflow:
	runner = inheritance.LineageWorkflow(
		dlimit = 0.03,
		flimit = 0.97,
		additive_cutoff = 0.03,
		subtractive_cutoff = 0.03,
		derivative_cutoff = 0.01
	)
	return runner

def test_periodic_selection(lineage_workflow):
	filename = DATA_FOLDER / "model_periodic_selection.tsv"
	table = dataio.import_table(filename, index = 'Genotype')

	expected_lineage = {
		'genotype-1': 'genotype-0',
		'genotype-2': 'genotype-1',
		'genotype-3': 'genotype-2',
		'genotype-4': 'genotype-3'
	}

	ancestry = lineage_workflow.run(table)

	logger.debug(ancestry.to_table().to_string())

	assert ancestry.as_dict() == expected_lineage

def test_clonal_interferance(lineage_workflow):
	pass
def test_strong_selection(lineage_workflow):
	pass
