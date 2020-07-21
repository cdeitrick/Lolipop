import argparse
import itertools
import math
from pathlib import Path
from typing import List, Optional, Union

try:
	from muller import dataio
except ModuleNotFoundError:
	from . import dataio

from dataclasses import dataclass

__VERSION__ = "0.8.1"
DEBUG = False


# For convienience. Helps with autocomplete.
# Not really used, will probably remove later.
@dataclass
class ProgramOptions(argparse.Namespace):
	filename: Path
	output_folder: Path
	sheetname: str = "Sheet1"
	fixed_breakpoint: Optional[float] = None
	detection_breakpoint: float = 0.03
	significant_breakpoint: float = 0.15
	mode: bool = False
	frequencies: List[float] = 0.10
	similarity_breakpoint: float = 0.05
	difference_breakpoint: float = 0.10
	is_genotype: bool = False
	use_filter: bool = True
	annotate_all: bool = False
	save_pvalue: bool = True
	use_strict_filter: bool = False
	method: str = 'matlab'
	metric: str = "similarity"
	known_genotypes: Optional[Path] = None


ACCEPTED_METHODS = ["matlab", "hierarchy", "twostep"]


# noinspection PyTypeChecker
def parse_workflow_options(program_options: argparse.Namespace) -> ProgramOptions:
	"""
		Generates values for each of the program parameters from the given parameters on the command line.
	Parameters
	----------
	program_options

	Returns
	-------

	"""
	#if program_options.name != 'lineage': return program_options
	if program_options.flimit is None:
		program_options.flimit = 1 - program_options.dlimit
	if program_options.known_genotypes:
		program_options.known_genotypes = Path(program_options.known_genotypes)
		starting_genotypes = dataio.import_file.parse_known_genotypes(program_options.known_genotypes)
	else:
		starting_genotypes = None

	program_options.starting_genotypes = starting_genotypes
	if program_options.output_folder is None:
		program_options.output_folder = Path(program_options.filename.parent)

	if program_options.similarity_cutoff is not None and program_options.similarity_cutoff < 0:
		program_options.similarity_cutoff = None

	return program_options


def _parse_frequency_option(frequency: Union[str, List[float]]) -> List[float]:
	if isinstance(frequency, str):
		if ',' in frequency:
			frequencies = list(map(float, frequency.split(',')))
		else:
			frequencies = float(frequency)
	elif isinstance(frequency, float):
		frequencies = frequency
	else:
		frequencies = [float(i) for i in frequency]
	if isinstance(frequencies, float):
		frequencies = [math.fsum(itertools.repeat(frequencies, i)) for i in range(int(1 / frequencies) + 1)]
		frequencies = [round(i, 2) for i in frequencies]
	if 0 not in frequencies:
		frequencies = [0] + frequencies
	frequencies = sorted(frequencies, reverse = True)
	return frequencies


#####################################################################################
################################ Custom Parsers ####################################
#####################################################################################

class FixedBreakpointParser(argparse.Action):
	def __call__(self, parser, namespace, values, option_string = None):
		detected = namespace.dlimit
		if values == "1":
			values = 1 - detected
		else:
			values = float(values)
		setattr(namespace, self.dest, values)


#####################################################################################
############################# Lineage Parser Groups #################################
#####################################################################################
def _create_parser_lineage_group_graphics(parser: argparse.ArgumentParser):
	##############################################################################################################################################
	# -------------------------------------------------------- Graphics Options ------------------------------------------------------------------
	##############################################################################################################################################

	group_graphics = parser.add_argument_group(title = "Graphics Options", description = "Options to refine how the output graphics are generated.")

	group_graphics.add_argument(
		"--annotate-all",
		help = "Adds all gene labels to the muller plots, instead of the top three.",
		action = "store_true",
		dest = "annotate_all"
	)
	group_graphics.add_argument(
		"--genotype-colors",
		help = "An optional map of genotypes to specified colors.",
		action = "store",
		type = Path,
		default = None,
		dest = "genotype_palette_filename"
	)
	group_graphics.add_argument(
		"--no-render",
		help = "The svg files will not be generated.",
		action = "store_false",
		dest = "render"
	)
	group_graphics.add_argument(
		"--no-outline",
		help = 'Disables the white outline in the muller plots.',
		action = 'store_false',
		dest = 'draw_outline'
	)
	group_graphics.add_argument(
		"--highlight",
		help = "A comma-separated list of genotype names or annotations to highlight in the generated graphics.",
		action = "store",
		dest = "highlight"
	)
	group_graphics.add_argument(
		"--highlight-color",
		help = "What the color of highlighted genotypes should be. Only HEX color codes are supported.",
		action = 'store',
		dest = 'highlight_color',
		default = "#FFFF00"
	)

	return group_graphics


# noinspection PyTypeChecker,PyTypeChecker,PyTypeChecker
def _create_parser_lineage_group_data(parser: argparse.ArgumentParser):
	##############################################################################################################################################
	# ----------------------------------------------------- Additional Input Files ---------------------------------------------------------------
	##############################################################################################################################################

	group_data = parser.add_argument_group(title = "Input Data", description = "Additional data that can be used to improve the analysis.")
	group_data.add_argument(
		"--filename-pairwise",
		help = "Path to a table with pairwise distance calculations from a previous run using identical input parameters. Should be located " \
			   "in `tables/.distances.tsv in the output folder generated from the previous run. This table will be used rather than re-calculating " \
			   "all the pairwise distances again which may take a long time for very large datasets.",
		action = "store",
		dest = "filename_pairwise",
		type = Path,
		default = None
	)

	group_data.add_argument(
		"--gene-alias",
		help = "An optional two-column file with more accurate gene names. This is usefull when using a reference annotated via prokka.",
		action = "store",
		type = Path,
		default = None,
		dest = 'alias_filename'
	)
	group_data.add_argument(
		"-g", "--known-genotypes",
		help = "A file with trajectories known to be in the same genotypes. "
			   "Each genotype is defined by a comma-delimited line with the labels of the member trajectories.",
		action = "store",
		default = None,
		dest = "known_genotypes"
	)
	group_data.add_argument(
		"--known-ancestry",
		help = "A file designating the known ancestry of certain genotypes. Formatted like the ggmuller edges table.",
		dest = 'known_ancestry',
		default = None,
		type = Path
	)


def _create_parser_lineage_group_genotype_generation(parser: argparse.ArgumentParser):
	##############################################################################################################################################
	# ------------------------------------------------------ General Analysis Options ------------------------------------------------------------
	##############################################################################################################################################
	analysis_group = parser.add_argument_group(title = "Genotype and Lineage Parameters")
	analysis_group.add_argument(
		"--threads",
		help = "The number of processes to use. Adding more threads than available cpu cores provides no speedup.",
		action = "store",
		dest = "threads",
		type = int,
		default = 1
	)

	analysis_group.add_argument(
		"--metric",
		help = "The distance metric to use when clustering mutaitons into genotypes.",
		action = "store",
		dest = "metric",
		type = str,
		default = "binomial"
	)
	analysis_group.add_argument(
		"--similarity-cutoff",
		help = "Used to group trajectories into genotypes. The 10th percentile of the trajectory distances is used if None.",
		action = "store",
		dest = "similarity_cutoff",
		type = float,
		default = None
	)
	analysis_group.add_argument(
		"-p", "--pvalue",
		help = "The p-value to use for the statistics tests",
		action = "store",
		dest = "pvalue",
		type = float,
		default = 0.05
	)

	analysis_group.add_argument(
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = FixedBreakpointParser,
		dest = 'flimit',
		type = float,
		default = None
	)
	analysis_group.add_argument(
		"-d", "--detection",
		help = "The uncertainty to apply when performing frequency-based calculations. \
			For example, a frequency at a given timepoint is considered undetected if it falls below 0 + `uncertainty`.",
		action = 'store',
		default = 0.03,
		dest = 'dlimit',
		type = float
	)
	analysis_group.add_argument(
		"-s", "--significant",
		help = "The frequency at which to consider a genotype significantly greater than zero.",
		action = 'store',
		default = 0.15,
		dest = "slimit",
		type = float
	)
	analysis_group.add_argument(
		"--covariance",
		help = "Controls how much a nested and unnested genotype should be correlated/anticorrelated to be considered significant",
		default = 0.01,
		dest = "derivative_cutoff",
		type = float
	)
	analysis_group.add_argument(
		"--liberal",
		help = "Whether to allow genotypes to be placed under more recent nested genotypes when an earlier genotype also passes the score check.",
		action = "store_false",
		dest = "conservative"
	)

	return analysis_group


def _create_parser_lineage_group_filter(parser: argparse.ArgumentParser):
	##############################################################################################################################################
	# --------------------------------------------------- Genotype Filtering Options -------------------------------------------------------------
	##############################################################################################################################################
	group_filter = parser.add_argument_group(title = "Filtering Parameters", description = "Parameters to control the filtering process.")
	group_filter.add_argument(
		"--disable-all-filters",
		help = "Disables genotype filtering.",
		action = 'store_false',
		dest = 'use_filter'
	)
	group_filter.add_argument(
		"--strict-filter",
		help = """By default, the filters allow trajectories to appear both before and after a genotype"""
			   """fixes as long as they were undetected at the timepoint the sweep occurs. This generally"""
			   """represents mutations which appear, are removed during a genotype sweep, and reappear """
			   """afterwards. Using `--strict-filter` would remove these trajectories.""",
		action = "store_true",
		dest = "use_strict_filter"
	)
	group_filter.add_argument(
		"--disable-filter-single",
		help = "Whether to remove trajectories which only exist at a single timepoint.",
		action = "store_false",
		dest = "use_filter_single"
	)
	group_filter.add_argument(
		"--disable-filter-startsfixed",
		help = "Whether to remove mutations which are 'fixed' at the first timepoint.",
		action = "store_false",
		dest = "use_filter_startsfixed"
	)
	group_filter.add_argument(
		"--filter-constant",
		help = "Sets the delta value by which a mutational trajectory must vary by to not be removed for being a constant mutation. Set to 0 to disable",
		action = "store",
		dest = "filter_constant",
		default = 0.1,
		type = float
	)
	group_filter.add_argument(
		"--disable-genotype-filter",
		help = "Disables filtering based on genotypes.",
		action = "store_false",
		dest = "use_filter_genotype"
	)

	return group_filter


# noinspection PyTypeChecker,PyTypeCheckerFalse
def _create_parser_lineage_group_main(parser: argparse.ArgumentParser):
	group_main = parser.add_argument_group(title = "Main Options")

	group_main.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename',
		type = Path,
		required = True
	)
	group_main.add_argument(
		'-o', '--output',
		help = "The folder to save the files to.",
		action = 'store',
		dest = 'output_folder',
		type = Path,
		default = None
	)
	group_main.add_argument(
		"--name",
		help = 'Prefix to use when naming the output files. defaults to the dataset filename.',
		action = 'store',
		type = str,
		dest = 'name',
		default = None
	)
	##############################################################################################################################################
	# --------------------------------------------------------- Input Data Options ---------------------------------------------------------------
	##############################################################################################################################################
	# input_group = parser.add_argument_group(title = "Input Dataset Options", description = "Options related to the input dataset.")

	group_main.add_argument(
		"--sheetname",
		help = "Indicates the sheet to use if the input table is an excel workbook and the data is not in Sheet1",
		action = 'store',
		dest = 'sheetname',
		default = 0
	)
	group_main.add_argument(
		"--genotypes", "--cohorts",
		help = "Indicates that the input table contains genotypes rather than trajectories.",
		action = 'store_true',
		dest = 'is_genotype'
	)


	return group_main


#####################################################################################
########################## Main Application Parsers #################################
#####################################################################################
def create_lineageplot_parser(subparsers):
	parser = subparsers.add_parser('lineageplot')
	group_main = parser.add_argument_group(title = "Lineageplot Options")
	group_main.add_argument(
		"edges",
		help = "Path to an `edges` table with both `identity` and `parent` columns. An optional `annotations` column"
			   "may also be included.",
		type = Path
	)
	group_main.add_argument(
		"filename",
		help = "Name of where to save the plot.",
		type = Path
	)
	group_main.add_argument(
		"--sheetname",
		help = "The sheet to use if the input file is an excel sheet.",
		type = str,
		required = False
	)
	return group_main


def create_lineage_parser(subparsers) -> argparse.ArgumentParser:
	parser = subparsers.add_parser("lineage")

	# Broke up the indivisual groups into their own functions because this function was large enough to make browsing it annoying.
	_create_parser_lineage_group_main(parser)
	_create_parser_lineage_group_data(parser)
	_create_parser_lineage_group_genotype_generation(parser)
	# _create_parser_group_filter(parser)
	_create_parser_lineage_group_graphics(parser)

	return parser


# noinspection PyTypeChecker,PyTypeChecker
def create_benchmark_parser(subparsers) -> argparse.ArgumentParser:
	""" Defines options for utilities that complement the scripts.
		TODO: Add a script to directly convert breseq output to muller input.
	"""
	parser_benchmark: argparse.ArgumentParser = subparsers.add_parser(
		"benchmark",
		help = "Runs a series of benchmarks to test the time required to calculate all pairwise distances of a dataset based on available processes."
	)
	parser_benchmark.add_argument(
		"--dataset",
		help = "The location of a test dataset to use. Should be large enough that adding additional processes will actually imrpove parsing speed.",
		type = Path
	)
	parser_benchmark.add_argument(
		"--output",
		help = "The filenme of the output graph. Should be a .png file.",
		type = Path
	)

	return parser_benchmark


# noinspection PyTypeChecker,PyTypeChecker,PyTypeChecker
def create_muller_parser(subparsers) -> argparse.ArgumentParser:
	parser_muller: argparse.ArgumentParser = subparsers.add_parser(
		"muller",
		help = "Generates muller plots from a pair of `population` and `edges` tables."
	)

	parser_muller.add_argument(
		"population",
		help = "Path to the population table.",
		type = Path
	)

	parser_muller.add_argument(
		"edges",
		help = "Path to the edges table delineating the lineage of the population genotypes.",
		type = Path
	)
	parser_muller.add_argument(
		"output",
		help = "Path to save the generated muller plot as.",
		type = Path
	)

	return parser_muller


def create_timeseries_parser(subparsers) -> argparse.ArgumentParser:
	parser_timeseries: argparse.ArgumentParser = subparsers.add_parser(
		"timeseries",
		help = "A small utility to plot timeseries tables."
	)
	parser_timeseries.add_argument(
		"timeseries",
		help = "A table with either a `Trajectory` column or a `Genotype` column.",
		type = Path
	)
	return parser_timeseries


def create_parser() -> argparse.ArgumentParser:
	parser_parent = argparse.ArgumentParser(
		description = "Infers genotypes and lineage based on mutation frequency timeseries data observed over the course of an evolutionary biology experiment.",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	parser_parent.add_argument(
		"-v", "--version",
		action = 'version',
		version = f"%(prog)s {__VERSION__}"
	)

	subparsers = parser_parent.add_subparsers(dest = 'name')  # Each subparser can be identifies by the `name` attribute.
	create_lineage_parser(subparsers)
	create_benchmark_parser(subparsers)
	create_muller_parser(subparsers)
	create_lineageplot_parser(subparsers)

	return parser_parent


def get_arguments(arguments: Optional[List[str]] = None) -> argparse.Namespace:
	""" Implemented here to make sure the default parameters are properly applied. """
	parser = create_parser()
	args = parser.parse_args(arguments)

	args = parse_workflow_options(args)
	return args

