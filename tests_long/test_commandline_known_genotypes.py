from pathlib import Path
from typing import *
import pytest
import subprocess
from muller.dataio import projectoutput
data_folder = Path(__file__).parent.parent / "tests_fixed_lineage"
from loguru import logger
def run_command(filename: Path, output_folder: Path, known_genotypes:Path):
	""" Runs lolipop using the table with path `filename`. """
	filename_lolipop_script = Path.home() / "Documents" / "github" / "muller_diagrams" / "lolipop"
	command = [
		"python", filename_lolipop_script,
		"lineage",
		"--input", filename,
		"--output", output_folder,
		"--known-genotypes", known_genotypes
	]
	command = [str(i) for i in command]
	logger.debug(" ".join(command))

	process = subprocess.run([str(i) for i in command], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	filename_stdout = output_folder / "stdout.txt"
	filename_stderr = output_folder / "stderr.txt"
	filename_stdout.write_bytes(process.stdout)
	filename_stderr.write_bytes(process.stderr)

	return projectoutput.OutputStructure(output_folder)


def test_known_genotypes_lang(tmp_path):

	known_genotypes = [
		["trajectory-aqua-1","trajectory-aqua-2", "trajectory-aqua-3", "trajectory-aqua-4", "trajectory-aqua-5", "trajectory-aqua-6"],
		["trajectory-green-1", "trajectory-green-2", "trajectory-green-3"],
		["trajectory-red-1"],
		["trajectory-gold-1", "trajectory-gold-2", "trajectory-gold-3"],
		["trajectory-sienna-1"],
		["trajectory-orange-1"]
	]
	filename_input_table = data_folder / "dataset_lang" / "input_table.tsv"
	output_folder = tmp_path / "lang_output"
	filename_known_genotypes =tmp_path/ "known_genotypes.txt"

	with filename_known_genotypes.open('w') as file1:
		for genotype in known_genotypes:
			line = ",".join(genotype)
			file1.write(f"{line}\n")

	result = run_command(filename_input_table, output_folder, filename_known_genotypes)

	assert len(result.table_genotypes) == len(known_genotypes)



