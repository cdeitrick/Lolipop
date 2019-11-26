from pathlib import Path
import pytest
from muller.dataio import generate_output
import json
import pandas
DATA_FOLDER = Path(__file__).parent.parent / "data"

@pytest.fixture
def output_folder()->Path:
	return DATA_FOLDER / "example2"

@pytest.fixture
def data_basic(output_folder)-> generate_output.DataWorkflow:

	program_options_filename = output_folder / "supplementary-files" / "traverse-etal-B1-mutationfrequencies.trajectories.options.json"

	result = generate_output.DataWorkflow(
		version = "0.7.0",
		filename = output_folder / "traverse-etal-B1-mutationfrequencies.xlsx",
		program_options = json.loads(program_options_filename.read_text()),
		output_folder = output_folder
	)
	return result

@pytest.fixture
def data_genotypes(output_folder)-> generate_output.DataGenotypeInference:
	filename_trajectories = output_folder / "tables" / "nature12344-s2.BYB1-G07.trajectories.original.tsv"
	filename_genotypes = output_folder / "nature12344-s2.BYB1-G07.genotypes.tsv"
	result = generate_output.DataGenotypeInference(
		info = None,
		original_trajectories = pandas.read_csv(filename_trajectories),
		table_genotypes = pandas.read_csv(filename_genotypes),
		genotype_members = None
	)

def test_works(data_basic):
	from pprint import pprint
	pprint(data_basic)

if __name__ == "__main__":
	pass