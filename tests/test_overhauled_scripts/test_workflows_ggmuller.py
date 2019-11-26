from pathlib import Path

from muller import dataio, workflows
import pytest
import pandas
from loguru import logger
@pytest.fixture
def workflow_gg(tmp_path)->workflows.GGMuller:
	return workflows.GGMuller(folder = tmp_path)


DATA_FOLDER = Path(__file__).parent.parent / "data"
ggmuller_tables = DATA_FOLDER / "tables_ggmuller"

def test_workflow_ggmuller_generate_filesystem_structure(tmp_path, workflow_gg):
	folder = tmp_path / "folder"

	folder_data = folder / "data"
	folder_figure = folder / "figures"
	expected_paths = {
		'outputTablePopulationLong': folder_data / "populations.long.tsv",
		'outputTablePopulationWide': folder_data / "populations.wide.tsv",
		'outputTableEdges':      folder_data / "edges.tsv",
		'outputTableMuller':     folder_data / "mullerformat.tsv",

		'outputFigureMuller':    folder_figure / "mullerplot.png"
	}

	result = workflow_gg.generate_filesystem_structure(folder)

	assert result == expected_paths

def test_workflow_ggmuller_import_table(workflow_gg):
	filename = ggmuller_tables / "B1_muller_try1.ggmuller.populations.tsv"

	result = workflow_gg.load_table(filename)
	# Testing without pandas.testing.assertframeequal because that method
	# does additional testing unrelated to what we want to use here (ex. testing df.index.name)

	assert len(result) == 84

	assert list(result.columns) == ["Generation","Identity", "Population"]

	assert sorted(set(result['Generation'].values)) == [0, 17,25,44,66,75,90]

def test_workflow_ggmuller_convert_to_timeseries(workflow_gg):
	pass


def test_workflow_ggmuller_with_default_paths(tmp_path):
	""" Tests whether the GGMuller workflow utilizes the initial paths"""
	# First, test the workflow using tables formatted for use by gmuller.
	filename_table_basic_populations = DATA_FOLDER / "tables_ggmuller" / "B1_muller_try1.ggmuller.populations.tsv"
	filename_table_basic_edges = DATA_FOLDER / "tables_ggmuller" / "B1_muller_try1.ggmuller.edges.tsv"

	output_folder = dataio.checkdir(tmp_path / "ggmuller_output")

	diagram_generator = workflows.GGMuller(folder = output_folder)

	# Test whether the `paths` variable is correct
	folder_data = output_folder / "data"
	folder_figure = output_folder / "figures"
	expected_paths = {
		'outputTablePopulationLong': folder_data / "populations.long.tsv",
		'outputTablePopulationWide': folder_data / "populations.wide.tsv",
		'outputTableEdges':      folder_data / "edges.tsv",
		'outputTableMuller':     folder_data / "mullerformat.tsv",
		'outputFigureMuller':    folder_figure / "mullerplot.png"
	}
	assert diagram_generator.paths == expected_paths

	result_paths = diagram_generator.run(filename_table_basic_populations, filename_table_basic_edges, output_folder)

	# Test whether each of those files were actually generated.
	for key, filename in expected_paths.items():
		logger.debug(f"Testing whether {key} exists...")
		assert filename.exists()
		logger.debug(f"'{key}' exists.")

	assert result_paths == expected_paths

	# Test whether the input tables were also saved.
	expected_table_edges = pandas.read_csv(filename_table_basic_edges, sep = '\t')

	result_table_edges = pandas.read_csv(result_paths['outputTableEdges'], sep = '\t')

	logger.debug(expected_table_edges.to_string())
	logger.debug(result_table_edges.to_string())

	expected_table_edges = expected_table_edges.set_index('Identity')['Parent']
	result_table_edges = result_table_edges.set_index('Identity')['Parent']
	pandas.testing.assert_series_equal(expected_table_edges, result_table_edges)
