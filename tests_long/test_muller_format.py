import itertools
import subprocess
from pathlib import Path

import pandas
import pytest
from loguru import logger

from muller.dataio.mullerformat import AdjacencyMatrix, GenerateMullerDataFrame


@pytest.fixture
def muller_dataframe_generator() -> GenerateMullerDataFrame:
	return GenerateMullerDataFrame()


FOLDER_DATA = Path(__file__).parent.parent / "tests" / "data"
FOLDER_TABLES = FOLDER_DATA / "tables_ggmuller"
FOLDER_SCRIPTS = FOLDER_DATA / "scripts_ggmuller"
tables_combined = [
	(FOLDER_TABLES / "B1_muller_try1.ggmuller.populations.tsv", FOLDER_TABLES / "B1_muller_try1.ggmuller.edges.tsv"),
	(FOLDER_TABLES / "B1_Muller.ggmuller.populations.tsv", FOLDER_TABLES / "B1_Muller.ggmuller.edges.tsv"),
	(FOLDER_TABLES / "Planktonic3.ggmuller.populations.tsv", FOLDER_TABLES / "Planktonic_3_mullerinput.ggmuller.edges.tsv"),
	(FOLDER_TABLES / "mullerinputbio1.ggmuller.populations.tsv", FOLDER_TABLES / "mullerinputbio1.ggmuller.edges.tsv"),
	(FOLDER_TABLES / "P1_Final_Muller.ggmuller.populations.tsv", FOLDER_TABLES / "P1_Final_Muller.ggmuller.edges.tsv"),
	# The LTEE tables are causing an issue where timepoints are getting duplicated.
	(FOLDER_TABLES / "m5_correct.ggmuller.populations.tsv", FOLDER_TABLES / "m5_correct.ggmuller.edges.tsv")
]

tables_population = [i[0] for i in tables_combined]


class ScriptFilenames:
	script_folder = FOLDER_SCRIPTS

	# Mainly used to collapse the filename methods in an IDE.
	@classmethod
	def get_script_full(cls) -> Path:
		return cls.script_folder / "rscript.adjustpopulations.r"

	@classmethod
	def get_script_find_start_positions(cls) -> Path:
		return cls.script_folder / "rscript.adjustpopulations.addstartpoints.findstartpositions.r"

	@classmethod
	def get_script_get_initial_generations(cls) -> Path:
		return cls.script_folder / "rscript.adjustpopulations.addstartpoints.getinitialgenerations.r"

	@classmethod
	def get_script_adjust_population(cls) -> Path:
		return cls.script_folder / "rscript.adjustpopulations.addstartpoints.adjustpopulation.r"

	@classmethod
	def get_script_run(cls) -> Path:
		return cls.script_folder / "rscript.adjustpopulations.addstartpoints.run.r"

	@classmethod
	def get_script_get_genotype_ages(cls) -> Path:
		return cls.script_folder / "rscript.getmullerdf.genotypeages.r"

	@classmethod
	def get_script_add_semi_frequencies(cls) -> Path:
		return cls.script_folder / "rscript.getmullerdf.addsemifrequencies.r"

	@classmethod
	def get_script_generate_unique_ids_genotypes(cls) -> Path:
		return cls.script_folder / "rscript.reorderbyvector.assigneuniqueids.r"

	@classmethod
	def get_script_generate_unique_ids_generations(cls) -> Path:
		return cls.script_folder / "rscript.reorderbyvector.generateuniqueidsgeneration.r"

	@classmethod
	def get_script_clean_edges(cls) -> Path:
		return cls.script_folder / "rscript.adjacencymatrix.cleanedges.r"

	@classmethod
	def get_script_moves(cls) -> Path:
		return cls.script_folder / "rscript.adjacencymatrix.move.r"

	@classmethod
	def get_script_find_start_node(cls) -> Path:
		return cls.script_folder / "rscript.adjacencymatrix.findstartnode.r"

	@classmethod
	def get_script_path_vector(cls) -> Path:
		return cls.script_folder / "rscript.adjacencymatrix.pathvector.r"

	@classmethod
	def get_script_generate_muller_dataframe(cls) -> Path:
		return cls.script_folder / "rscript.generatemullerdf.r"


def approxlists(left, right) -> bool:
	# Used to test whether two lists of floats are identical. Using `==` will not work with floats sometimes due to precision.
	result = [pytest.approx(i, j) for i, j in zip(left, right)]
	return all(result)


def get_lookup_table(population_filename) -> pandas.Series:
	abs_path = FOLDER_DATA / population_filename
	mullerdf = GenerateMullerDataFrame()

	population = pandas.read_csv(abs_path, sep = '\t')
	population_converted = mullerdf.apply_add_start_points(population)
	population_generations = mullerdf.get_genotype_ages(population_converted)

	return population_generations


def get_muller_df(filename: Path) -> pandas.DataFrame:
	df = pandas.read_csv(filename, sep = "\t")
	return GenerateMullerDataFrame().correct_population_values(df)


def checkdir(path: Path) -> Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path


@pytest.mark.parametrize("populationfilename", tables_population)
def test_add_start_points(tmp_path, muller_dataframe_generator, populationfilename):
	filename_truthset = tmp_path / "truth"
	population = pandas.read_csv(populationfilename, sep = '\t')
	command = ["Rscript", ScriptFilenames.get_script_full(), populationfilename, filename_truthset]
	subprocess.run(command)
	truthset = pandas.read_csv(filename_truthset, sep = "\t")
	truthset['Generation'] = truthset['Generation'].astype(float)

	result = muller_dataframe_generator.apply_add_start_points(population)

	# pandas.testing does other stuff we dont care about.

	pandas.testing.assert_frame_equal(
		result.reset_index(), truthset.reset_index(),
		check_exact = False,
		check_less_precise = True
	)


@pytest.mark.parametrize("filename_population", tables_population)
def test_find_start_points(tmp_path, muller_dataframe_generator, filename_population: Path):
	filename_truthset = tmp_path / "truthset"
	population = pandas.read_csv(filename_population, sep = '\t')

	command = ["Rscript", ScriptFilenames.get_script_find_start_positions(), filename_population, filename_truthset]
	subprocess.run(command)

	truthset = filename_truthset.read_text()
	truthset = float(truthset)
	result = muller_dataframe_generator._find_start_positions(population['Generation'], 0.5)
	assert result == truthset


@pytest.mark.parametrize("filename_population", tables_population)
def test_get_initial_generations(tmp_path, muller_dataframe_generator, filename_population):
	filename_truthset = tmp_path / "truth"
	population = pandas.read_csv(filename_population, sep = '\t')
	command = ["Rscript", ScriptFilenames.get_script_get_initial_generations(), filename_population, filename_truthset]
	subprocess.run(command)

	truthset = pandas.read_csv(filename_truthset, sep = '\t')

	result = muller_dataframe_generator._get_initial_generations(population)

	pandas.testing.assert_frame_equal(result.reset_index(drop = True), truthset.reset_index(drop = True))


@pytest.mark.parametrize("filename_population", tables_population)
def test_adjust_population_table(tmp_path, muller_dataframe_generator, filename_population):
	filename_truth = tmp_path / "truthset"
	population = pandas.read_csv(filename_population, sep = '\t')

	command = ["Rscript", ScriptFilenames.get_script_adjust_population(), filename_population, filename_truth]
	subprocess.run(command)
	truthset = pandas.read_csv(filename_truth, sep = '\t')
	truthset['Generation'] = truthset['Generation'].astype(float)
	# Make sure the truthset generations column is `float`
	truthset['Generation'] = truthset['Generation'].astype(float)

	first_gens = muller_dataframe_generator._get_initial_generations(population)
	result = muller_dataframe_generator._adjust_population_table(population, first_gens, 0.5)

	pandas.testing.assert_frame_equal(result.reset_index(drop = True), truthset.reset_index(drop = True))


@pytest.mark.parametrize("filename_population", tables_population)
def test_add_start_points_run(tmp_path, muller_dataframe_generator, filename_population):
	filename_truthset = tmp_path / "truthset"
	population = pandas.read_csv(filename_population, sep = '\t')
	command = ["Rscript", ScriptFilenames.get_script_run(), filename_population, filename_truthset]
	subprocess.run(command)
	truthset = pandas.read_csv(filename_truthset, sep = '\t')
	truthset['Generation'] = truthset['Generation'].astype(float)  # Will be int sometimes.

	result = muller_dataframe_generator.apply_add_start_points(population)

	pandas.testing.assert_frame_equal(result, truthset)


@pytest.mark.parametrize("filename_population", tables_population)
def test_get_genotype_ages(tmp_path, muller_dataframe_generator, filename_population):
	filename_sorted = tmp_path / "pop_sorted"
	filename_truth = tmp_path / "truthset"

	population = pandas.read_csv(filename_population, sep = '\t')
	# Simpler for now just to rerun the process, since it's already been checked above.
	pop_df_adjusted = muller_dataframe_generator.apply_add_start_points(population)
	pop_df_sorted = pop_df_adjusted.sort_values(by = ["Generation", "Population", "Identity"], ascending = [True, False, True])
	pop_df_sorted.to_csv(filename_sorted, sep = "\t", index = False)
	result_ages = muller_dataframe_generator.get_genotype_ages(pop_df_sorted).reset_index()

	command = ["Rscript", ScriptFilenames.get_script_get_genotype_ages(), filename_sorted, filename_truth]
	subprocess.run(command)
	truthset = pandas.read_csv(filename_truth, sep = "\t")

	pandas.testing.assert_frame_equal(result_ages.reset_index(drop = True), truthset.reset_index(drop = True))


@pytest.mark.parametrize("population_filename", tables_population)
def test_add_semi_frequencies(tmp_path, muller_dataframe_generator, population_filename):
	filename_truthset = tmp_path / "truth"
	filename_converted = tmp_path / "population_converted"
	population = pandas.read_csv(population_filename, sep = '\t')
	population_converted = muller_dataframe_generator.apply_add_start_points(population)

	population_converted.to_csv(filename_converted, sep = "\t", index = False)
	population_rank = muller_dataframe_generator.add_semi_frequencies(population_converted)

	command = ["Rscript", ScriptFilenames.get_script_add_semi_frequencies(), filename_converted, filename_truthset]
	subprocess.run(command)
	truthset = pandas.read_csv(filename_truthset, sep = '\t')
	truthset['Generation'] = truthset['Generation'].astype(float)
	pandas.testing.assert_frame_equal(population_rank.reset_index(drop = True), truthset.reset_index(drop = True))


@pytest.mark.parametrize("filename_population,filename_edges", tables_combined)
def test_clean_edges(filename_population, filename_edges, tmp_path):
	filename_lookup = tmp_path / "lookup"
	filename_truthset = tmp_path / "truthset"
	logger.debug(filename_truthset)

	lookup_table = get_lookup_table(filename_population)
	lookup_table.to_frame().to_csv(filename_lookup, sep = "\t")
	table_edges = pandas.read_csv(filename_edges, sep = "\t")

	edges = AdjacencyMatrix(table_edges)
	result = edges.generate_clean_edges(lookup_table)

	command = ["Rscript", ScriptFilenames.get_script_clean_edges(), filename_edges, filename_lookup, filename_truthset]
	subprocess.run(command)
	truthset = pandas.read_csv(filename_truthset, sep = "\t")

	pandas.testing.assert_frame_equal(result.reset_index(drop = True), truthset.reset_index(drop = True))


@pytest.mark.parametrize("filename_population,filename_edges", tables_combined)
def test_moves(tmp_path, filename_population, filename_edges):
	filename_am = tmp_path / "adjacencymatrix"

	table_edges = pandas.read_csv(filename_edges, sep = "\t")
	table_lookup = get_lookup_table(filename_population)

	am = AdjacencyMatrix(table_edges, table_lookup)
	am.clean_edges.to_csv(filename_am, sep = "\t", index = False)

	command = ["Rscript", ScriptFilenames.get_script_moves(), filename_am]

	# Test move_up

	for index in am.clean_edges['Identity'].values:
		output_filename = tmp_path / f"output.{filename_edges.stem}.{index}"
		move_command = command + [str(index), output_filename]
		logger.debug(move_command)
		subprocess.run(move_command)
		truth = output_filename.read_text()
		truth = [float(i) for i in truth.split(' ')]

		moveup = am.move_up(index)
		movedown = am.move_down(index)
		moveright = am.move_right(index)

		result = [moveup, movedown, moveright]

		assert result == truth


@pytest.mark.parametrize("filename_population, filename_edges", tables_combined)
def test_find_start_node(tmp_path, filename_population, filename_edges):
	filename_clean_edges = tmp_path / "clean_edges"
	filename_truth_node = tmp_path / "truth_node"

	table_edges = pandas.read_csv(filename_edges, sep = "\t")
	table_lookup = get_lookup_table(filename_population)

	am = AdjacencyMatrix(table_edges, table_lookup)
	result = am.find_start_node()
	am.clean_edges.to_csv(filename_clean_edges, sep = "\t", index = False)

	command = ["Rscript", ScriptFilenames.get_script_find_start_node(), filename_clean_edges, filename_truth_node]

	subprocess.run(command)

	truth_value = filename_truth_node.read_text()
	truth_value = float(truth_value)

	assert result == truth_value


@pytest.mark.parametrize("filename_population, filename_edges", tables_combined)
def test_path_vector(tmp_path, filename_population, filename_edges):
	filename_clean_edges = tmp_path / "clean_edges"
	filename_truth_result = tmp_path / "truthset"

	table_edges = pandas.read_csv(filename_edges, sep = "\t")
	table_lookup = get_lookup_table(filename_population)

	am = AdjacencyMatrix(table_edges, table_lookup)
	am.clean_edges.to_csv(filename_clean_edges, sep = "\t", index = False)
	result = am.path_vector()
	result = [am.ages.loc[i] for i in result]  # Need to convert back to ages since that is what the rscripts expect.

	command = ["Rscript", ScriptFilenames.get_script_path_vector(), filename_clean_edges, filename_truth_result]
	subprocess.run(command)
	# Extract the truth value. R automatically inserts newlines, so have to filter those out.
	truth_text = filename_truth_result.read_text()
	lines = [i.strip().split(" ") for i in truth_text.split("\n")]
	truth_value = list(int(i) for i in itertools.chain.from_iterable(lines) if i)

	assert result == truth_value


def test_reorder_by_vector_generate_unique_ids_genotype(tmp_path):
	filename_script = ScriptFilenames().get_script_generate_unique_ids_genotypes()
	filename_truth = tmp_path / "truthvalue"
	vector = [f"genotype-{i}" for i in [0, 5, 3, 3, 2, 6, 6, 2, 1, 7, 7, 1, 5, 4, 4, 0]]
	logger.info(vector)
	result = GenerateMullerDataFrame()._reorder_by_vector_generate_unique_ids_genotype(vector)

	command = ["Rscript", filename_script, ",".join(vector), filename_truth]
	subprocess.run(command)

	truth_table = pandas.read_csv(filename_truth, sep = "\t")
	assert result == truth_table["Unique_id"].tolist()


@pytest.mark.parametrize("filename_population, filename_edges", tables_combined[:-1])  # Known bug in the source ggmuller scripts.
def test_generate_muller_dataframe(tmp_path, filename_population, filename_edges):
	""" Tests whether the full get_Muller_df function has been implemented correctly in python."""
	filename_script = ScriptFilenames.get_script_generate_muller_dataframe()
	filename_truthset = tmp_path / "truthset"

	muller_generator = GenerateMullerDataFrame()

	table_edges = pandas.read_csv(filename_edges, sep = "\t")
	table_population = pandas.read_csv(filename_population, sep = "\t")
	result = muller_generator.run(table_edges, table_population)

	command = ["Rscript", filename_script, filename_edges, filename_population, filename_truthset]
	subprocess.run(command)

	truthset = pandas.read_csv(filename_truthset, sep = "\t")
	logger.debug(result.columns)
	logger.debug(truthset.columns)
	pandas.testing.assert_frame_equal(result.reset_index(drop = True), truthset.reset_index(drop = True))
