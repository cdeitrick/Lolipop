import argparse
import itertools
import math
import sys
from pathlib import Path
from typing import List, Optional, Union

try:
	from muller import dataio
except ModuleNotFoundError:
	from . import dataio

from dataclasses import dataclass, fields

__VERSION__ = "0.6.0"
DEBUG = False


# For convienience. Helps with autocomplete.
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

	def show(self):
		for field in fields(self):
			print(field)


ACCEPTED_METHODS = ["matlab", "hierarchy", "twostep"]


# Custom method added to argparse.ArgumentParser
def set_default_subparser(self, name, args = None, positional_args = 0):
	"""default subparser selection. Call after setup, just before parse_args()
	name: is the name of the subparser to call by default
	args: if set is the argument list handed to parse_args()

	, tested with 2.7, 3.2, 3.3, 3.4
	it works with 2.6 assuming argparse is installed
	"""

	subparser_found = False
	for arg in sys.argv[1:]:
		if arg in ['-h', '--help']:  # global help if no subparser
			break
	else:
		for x in self._subparsers._actions:
			if not isinstance(x, argparse._SubParsersAction):
				continue
			for sp_name in x._name_parser_map.keys():
				if sp_name in sys.argv[1:]:
					subparser_found = True
		if not subparser_found:
			# insert default in last position before global positional
			# arguments, this implies no global options are specified after
			# first positional argument
			if args is None:
				sys.argv.insert(len(sys.argv) - positional_args, name)
			else:
				args.insert(len(args) - positional_args, name)


def parse_workflow_options(program_options: ProgramOptions) -> ProgramOptions:
	"""
		Generates values for each of the program parameters from the given parameters on the command line.
	Parameters
	----------
	program_options

	Returns
	-------

	"""

	if program_options.fixed_breakpoint is None:
		program_options.fixed_breakpoint = 1 - program_options.detection_breakpoint
	if program_options.additive_cutoff is None:
		program_options.additive_cutoff = program_options.detection_breakpoint
	if program_options.subtractive_cutoff is None:
		program_options.subtractive_cutoff = program_options.detection_breakpoint

	if program_options.known_genotypes:
		program_options.known_genotypes = Path(program_options.known_genotypes)
		starting_genotypes = dataio.parse_known_genotypes(program_options.known_genotypes)
	else:
		starting_genotypes = None

	program_options.starting_genotypes = starting_genotypes
	cluster_method = program_options.method
	if cluster_method not in ACCEPTED_METHODS:
		message = f"{cluster_method} is not a valid option for the --method option. Expected one of {ACCEPTED_METHODS}"
		raise ValueError(message)
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


class FrequencyParser(argparse.Action):
	def __call__(self, parser, namespace, values, option_string = None):
		if not values:
			parsed_frequencies = self.default
		else:
			parsed_frequencies = _parse_frequency_option(values)
		setattr(namespace, self.dest, parsed_frequencies)


class FixedBreakpointParser(argparse.Action):
	def __call__(self, parser, namespace, values, option_string = None):
		detected = namespace.detection_breakpoint
		if values == "1":
			values = 1 - detected
		else:
			values = float(values)
		setattr(namespace, self.dest, values)


# noinspection PyTypeChecker

def create_utility_parser(subparsers) -> argparse.ArgumentParser:
	""" Defines options for utilities that complement the scripts.
		TODO: Add a script to directly convert breseq output to muller input.
	"""
	parser_benchmark: argparse.ArgumentParser = subparsers.add_parser(
		"benchmark",
		help = "Runs a series of benchmarks to test the time required to calculate all pairwise distances of a dataset based on available processes."
	)
	parser_benchmark.add_argument(
		"--benchmark",
		help = "",
		action = "store_true",
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


def _create_parser_group_graphics(parser: argparse.ArgumentParser):
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


def _create_parser_group_data(parser: argparse.ArgumentParser):
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


def _create_parser_group_analysis(parser: argparse.ArgumentParser):
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
		default = 2
	)
	analysis_group.add_argument(
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = FixedBreakpointParser,
		dest = 'fixed_breakpoint',
		type = float
	)
	analysis_group.add_argument(
		"-d", "--detection",
		help = "The uncertainty to apply when performing frequency-based calculations. \
			For example, a frequency at a given timepoint is considered undetected if it falls below 0 + `uncertainty`.",
		action = 'store',
		default = 0.03,
		dest = 'detection_breakpoint',
		type = float
	)
	analysis_group.add_argument(
		"-s", "--significant",
		help = "The frequency at which to consider a genotype significantly greater than zero.",
		action = 'store',
		default = 0.15,
		dest = "significant_breakpoint",
		type = float
	)
	analysis_group.add_argument(
		"--additive",
		help = "Controls how the additive score between a nested and unnested genotype is calculated. Defaults to the detection cutoff value.",
		action = 'store',
		default = None,
		dest = 'additive_cutoff',
		type = float
	)
	analysis_group.add_argument(
		"--subtractive",
		help = "Controls when the combined frequencies of a nested and unnested genotype are considered consistently larger than the fixed cutoff."
			   "Defaults to the detection cutoff value.",
		default = None,
		dest = "subtractive_cutoff"
	)
	analysis_group.add_argument(
		"--covariance",
		help = "Controls how much a nested and unnested genotype should be correlated/anticorrelated to be considered significant",
		default = 0.01,
		dest = "derivative_cutoff",
		type = float
	)

	##############################################################################################################################################
	# ----------------------------------------------- Options for individual analysis steps ------------------------------------------------------
	##############################################################################################################################################

	analysis_group.add_argument(
		"-f", "--frequencies",
		help = 'The frequency cutoff to use when sorting the muller_genotypes by first detected frequency. For example, a value of 0.15 will use the frequencies 0,.15,.30,.45...',
		action = FrequencyParser,
		dest = 'frequencies',
		default = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
	)
	analysis_group.add_argument(
		"--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint",
		type = float
	)
	analysis_group.add_argument(
		"--difference-cutoff",
		help = "Minimum p-value to consider a pair of genotypes unrelated. Used when splitting muller_genotypes when using the two-step method.",
		action = "store",
		default = 0.25,
		dest = "difference_breakpoint",
		type = float
	)
	return analysis_group


def _create_parser_group_filter(parser: argparse.ArgumentParser):
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

	return group_filter


def _create_parser_group_main(parser: argparse.ArgumentParser):
	group_main = parser.add_argument_group(title = "Main Options")

	group_main.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename',
		type = Path,
		# required = True
	)
	group_main.add_argument(
		'-o', '--output',
		help = "The folder to save the files to.",
		action = 'store',
		dest = 'output_folder',
		type = Path,
		# required = True
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


def _create_parser_group_clustering(parser: argparse.ArgumentParser):
	##############################################################################################################################################
	# --------------------------------------------------- Genotype Clustering Options ------------------------------------------------------------
	##############################################################################################################################################
	group_cluster = parser.add_argument_group(title = "Genotype Clustering Parameters", description = "Configures how genotypes are clustered.")
	group_cluster.add_argument(
		'-m', '--method',
		help = "The clustering method to use. `matlab` will use the original two-step algorithm while `hierarchy` will use hierarchical clustering.",
		action = "store",
		default = "hierarchy",
		dest = "method",
		choices = ['matlab', 'hierarchy', 'twostep']
	)
	group_cluster.add_argument(
		"--metric",
		help = "Selects the distance metric to use. Each metric tends to focus on a specific feature between two series, " \
			   "such as the difference between them or how well they are correlated.",
		action = "store",
		default = "binomial",
		dest = "metric",
		choices = ['similarity', 'binomial', 'pearson', 'minkowski', 'jaccard', 'combined']
	)


def create_parser() -> argparse.ArgumentParser:
	# Add the custom method to select the default parser
	argparse.ArgumentParser.set_default_subparser = set_default_subparser

	parser_parent = argparse.ArgumentParser(
		description = "Infers genotypes and lineage based on mutation frequency timeseries data observed over the course of an evolutionary biology experiment.",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)

	subparsers = parser_parent.add_subparsers()
	parser = subparsers.add_parser("lineage")
	parser.add_argument(
		"-v", "--version",
		action = 'version',
		version = f"%(prog)s {__VERSION__}"
	)

	# Broke up the indivisual groups into their own functions because this function was large enough to make browsing it annoying.
	_create_parser_group_main(parser)
	_create_parser_group_clustering(parser)
	_create_parser_group_data(parser)
	_create_parser_group_analysis(parser)
	_create_parser_group_filter(parser)
	_create_parser_group_graphics(parser)
	create_utility_parser(subparsers)
	return parser_parent
