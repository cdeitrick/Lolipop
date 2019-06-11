import argparse
import itertools
import math
from pathlib import Path
from typing import List, Optional, Union

try:
	from muller import dataio
except ModuleNotFoundError:
	from . import dataio

from dataclasses import dataclass, fields

__VERSION__ = "0.5.1"
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


def create_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser(
		description = "Generates muller diagrams based on a set of mutational trajectories.",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter
	)
	##############################################################################################################################################
	# --------------------------------------------------------- Required Parameters --------------------------------------------------------------
	##############################################################################################################################################
	parser.add_argument(
		"-v", "--version",
		action = 'version',
		version = f"%(prog)s {__VERSION__}"
	)
	parser.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename',
		type = Path,
		required = True
	)
	parser.add_argument(
		'-o', '--output',
		help = "The folder to save the files to.",
		action = 'store',
		dest = 'output_folder',
		type = Path,
		required = True
	)
	parser.add_argument(
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
	parser.add_argument(
		"--sheetname",
		help = "Indicates the sheet to use if the input table is an excel workbook and the data is not in Sheet1",
		action = 'store',
		dest = 'sheetname',
		default = 0
	)
	parser.add_argument(
		"--genotypes", "--cohorts",
		help = "Indicates that the input table contains genotypes rather than trajectories.",
		action = 'store_true',
		dest = 'is_genotype'
	)

	##############################################################################################################################################
	# ------------------------------------------------------ General Analysis Options ------------------------------------------------------------
	##############################################################################################################################################
	parser.add_argument(
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = FixedBreakpointParser,
		dest = 'fixed_breakpoint',
		type = float
	)
	parser.add_argument(
		"-d", "--detection",
		help = "The uncertainty to apply when performing frequency-based calculations. \
			For example, a frequency at a given timepoint is considered undetected if it falls below 0 + `uncertainty`.",
		action = 'store',
		default = 0.03,
		dest = 'detection_breakpoint',
		type = float
	)
	parser.add_argument(
		"-s", "--significant",
		help = "The frequency at which to consider a genotype significantly greater than zero.",
		action = 'store',
		default = 0.15,
		dest = "significant_breakpoint",
		type = float
	)
	parser.add_argument(
		"--additive",
		help = "Controls how the additive score between a nested and unnested genotype is calculated. Defaults to the detection cutoff value.",
		action = 'store',
		default = None,
		dest = 'additive_cutoff',
		type = float
	)
	parser.add_argument(
		"--subtractive",
		help = "Controls when the combined frequencies of a nested and unnested genotype are considered consistently larger than the fixed cutoff."
			   "Defaults to the detection cutoff value.",
		default = None,
		dest = "subtractive_cutoff"
	)
	parser.add_argument(
		"--derivative",
		help = "Controls how much a nested and unnested genotype should be correlated/anticorrelated to be considered significant",
		default = 0.01,
		dest = "derivative_cutoff",
		type = float
	)

	##############################################################################################################################################
	# ----------------------------------------------- Options for individual analysis steps ------------------------------------------------------
	##############################################################################################################################################

	parser.add_argument(
		"-f", "--frequencies",
		help = 'The frequency cutoff to use when sorting the muller_genotypes by first detected frequency. For example, a value of 0.15 will use the frequencies 0,.15,.30,.45...',
		action = FrequencyParser,
		dest = 'frequencies',
		default = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
	)
	parser.add_argument(
		"--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint",
		type = float
	)
	parser.add_argument(
		"--difference-cutoff",
		help = "Minimum p-value to consider a pair of muller_genotypes unrelated. Used when splitting muller_genotypes.",
		action = "store",
		default = 0.25,
		dest = "difference_breakpoint",
		type = float
	)

	##############################################################################################################################################
	# --------------------------------------------------- Genotype Clustering Options ------------------------------------------------------------
	##############################################################################################################################################
	parser.add_argument(
		"--no-filter",
		help = "Disables genotype filtering.",
		action = 'store_false',
		dest = 'use_filter'
	)
	parser.add_argument(
		"--include-single",
		help = "Whether to disable the filter which excludes trajectories which only exist at a single timepoint.",
		action = "store_false",
		dest = "include_single"
	)
	parser.add_argument(
		"--strict-filter",
		help = """By default, the filters allow trajectories to appear both before and after a genotype"""
			   """fixes as long as they were undetected at the timepoint the sweep occurs. This generally"""
			   """represents mutations which appear, are removed during a genotype sweep, and reappear """
			   """afterwards. Using `--strict-filter` would remove these trajectories.""",
		action = "store_true",
		dest = "use_strict_filter"
	)
	parser.add_argument(
		'-m', '--method',
		help = "The clustering method to use. `matlab` will use the original two-step algorithm while `hierarchy` will use hierarchical clustering.",
		action = "store",
		default = "hierarchy",
		dest = "method",
		choices = ['matlab', 'hierarchy', 'twostep']
	)
	parser.add_argument(
		"--metric",
		help = "Selects the distance metric to use. Each metric tends to focus on a specific feature between two series, " \
			   "such as the difference between them or how well they are correlated.",
		action = "store",
		default = "binomial",
		dest = "metric",
		choices = ['similarity', 'binomial', 'pearson', 'minkowski', 'jaccard', 'combined']
	)
	##############################################################################################################################################
	# -------------------------------------------------------- Graphics Options ------------------------------------------------------------------
	##############################################################################################################################################
	parser.add_argument(
		"--annotate-all",
		help = "Adds all gene labels to the muller plots, instead of the top three.",
		action = "store_true",
		dest = "annotate_all"
	)
	parser.add_argument(
		"--genotype-colors",
		help = "An optional map of genotypes to specified colors.",
		action = "store",
		type = Path,
		default = None,
		dest = "genotype_palette_filename"
	)

	##############################################################################################################################################
	# ----------------------------------------------------- Additional Input Files ---------------------------------------------------------------
	##############################################################################################################################################
	parser.add_argument(
		"--gene-alias",
		help = "An optional two-column file with more accurate gene names. This is usefull when using a reference annotated via prokka.",
		action = "store",
		type = Path,
		default = None,
		dest = 'alias_filename'
	)
	parser.add_argument(
		"-g", "--known-genotypes",
		help = "A file with trajectories known to be in the same genotypes. "
			   "Each genotype is defined by a comma-delimited line with the labels of the member trajectories.",
		action = "store",
		default = None,
		dest = "known_genotypes"
	)
	parser.add_argument(
		"--known-ancestry",
		help = "A file designating the known ancestry of certain genotypes. Formatted like the ggmuller edges table.",
		dest = 'known_ancestry',
		default = None,
		type = Path
	)


	##############################################################################################################################################
	# ------------------------------------------------------- Graphics Options -------------------------------------------------------------------
	##############################################################################################################################################

	parser.add_argument(
		"--no-outline",
		help = 'Disables the white outline in the muller plots.',
		action = 'store_true'
	)
	return parser
