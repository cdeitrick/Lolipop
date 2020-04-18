"""
	Runs the lineage workflow on the traverse dataset and compares it against the current lineages.

"""

from pathlib import Path
from typing import *
import subprocess
import pytest
import pandas
from loguru import logger
from muller.dataio import projectpaths
data_folder = Path(__file__).parent



class ExpectedStructure:
	def __init__(self, folder: Path):
		self.filename_genotypes = folder / "expected.genotypes.tsv"
		self.filename_edges = folder / "expected.edges.tsv"

		self.table_genotypes = pandas.read_csv(self.filename_genotypes, sep = "\t")
		self.table_edges = pandas.read_csv(self.filename_edges, sep = "\t")


def run_command(filename: Path, output_folder: Path):
	""" Runs lolipop using the table with path `filename`. """
	filename_lolipop_script = Path.home() / "Documents" / "github" / "muller_diagrams" / "lolipop"
	command = [
		"python", filename_lolipop_script,
		"lineage",
		"--input", filename,
		"--output", output_folder
	]

	process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	filename_stdout = output_folder / "stdout.txt"
	filename_stderr = output_folder / "stderr.txt"
	filename_stdout.write_bytes(process.stdout)
	filename_stderr.write_bytes(process.stderr)

	return projectpaths.OutputFilenames(output_folder, "input_table")


def dataset_test(folder_dataset: Path, output_folder: Path):
	filename_input_table = folder_dataset / "input_table.tsv"

	run_command(filename_input_table, output_folder)


def test_traverse_lineage(tmp_path):
	folder_traverse = data_folder / "dataset_traverse"
	folder_output = tmp_path / "test_dataset_traverse"

	expected = ExpectedStructure(folder_traverse)
	filename_input_table = folder_traverse / "input_table.tsv"

	output = run_command(filename_input_table, folder_output)

	# Test the genotype table
	output_table = pandas.read_csv(output.filename_table_genotypes, sep = "\t")
	pandas.testing.assert_frame_equal(expected.table_genotypes, output_table)

	# Test the edges table
	output_edges = pandas.read_csv(output.filename_table_edges, sep = "\t")
	pandas.testing.assert_frame_equal(expected.table_edges, output_edges)


def test_lang_lineage(tmp_path):
	folder_lang = data_folder / "dataset_lang"
	folder_output = tmp_path / "test_dataset_lang"
	filename_input_table = folder_lang / "input_table.tsv"

	expected = ExpectedStructure(folder_lang)

	output = run_command(filename_input_table, folder_output)

	# Test the genotype table
	pandas.testing.assert_frame_equal(expected.table_genotypes, output.table_genotypes)

	# Test the edges table
	pandas.testing.assert_frame_equal(expected.table_genotypes, output.table_genotypes)
