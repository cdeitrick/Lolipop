import argparse
import itertools
import math
from pathlib import Path
from typing import List, Optional, Union, Dict
import csv
try:
	from muller.options import GenotypeOptions, SortOptions, OrderClusterParameters
except ModuleNotFoundError:
	from options import GenotypeOptions, SortOptions, OrderClusterParameters

from dataclasses import dataclass, fields


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


ACCEPTED_METHODS = ["matlab", "hierarchy"]


def parse_workflow_options(program_options: ProgramOptions):
	# program_options = ProgramOptions.from_parser(program_options)
	if program_options.fixed_breakpoint is None:
		program_options.fixed_breakpoint = 1 - program_options.detection_breakpoint
	compatibility_mode = program_options.mode
	if compatibility_mode:
		program_options_genotype = GenotypeOptions.from_matlab()
		program_options_sort = SortOptions.from_matlab()
		program_options_clustering = OrderClusterParameters.from_matlab()
	else:
		if program_options.known_genotypes:
			program_options.known_genotypes = Path(program_options.known_genotypes)
			starting_genotypes = program_options.known_genotypes.read_text().split('\n')
			starting_genotypes = [i.split(',') for i in starting_genotypes if i]
		else:
			starting_genotypes = None
		program_options_genotype = GenotypeOptions(
			detection_breakpoint = program_options.detection_breakpoint,
			fixed_breakpoint = program_options.fixed_breakpoint,
			similarity_breakpoint = program_options.similarity_breakpoint,
			difference_breakpoint = program_options.difference_breakpoint,
			n_binom = None,
			method = program_options.method,
			metric = program_options.metric if program_options.method != 'similarity' else 'binomialp',
			starting_genotypes = starting_genotypes
		)
		program_options_clustering = OrderClusterParameters.from_breakpoints(
			program_options.detection_breakpoint,
			program_options.significant_breakpoint
		)

		program_options_sort = SortOptions(
			detection_breakpoint = program_options_genotype.detection_breakpoint,
			fixed_breakpoint = program_options_genotype.fixed_breakpoint,
			significant_breakpoint = program_options.significant_breakpoint,
			frequency_breakpoints = program_options.frequencies
		)
	cluster_method = program_options.method
	if cluster_method not in ACCEPTED_METHODS:
		message = f"{cluster_method} is not a valid option for the --method option. Expected one of {ACCEPTED_METHODS}"
		raise ValueError(message)
	return program_options, program_options_genotype, program_options_sort, program_options_clustering


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
		print("PARSING FIXED BREAKPOINT", values)
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
	parser.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename',
		type = Path,
		required = True
	)
	parser.add_argument(
		"--sheetname",
		help = "Indicates the sheet to use if the input table is an excel workbook and the data is not in Sheet1",
		action = 'store',
		dest = 'sheetname',
		default = 'Sheet1',
		type = str
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
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = FixedBreakpointParser,
		dest = 'fixed_breakpoint',
		type = float
	)
	parser.add_argument(
		"-u", "--uncertainty",
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
		"--matlab",
		help = "Mimics the output of the original matlab script.",
		action = 'store_true',
		dest = "mode"
	)
	parser.add_argument(
		"-f", "--frequencies",
		help = 'The frequency cutoff to use when sorting the muller_genotypes by first detected frequency. For example, a value of 0.15 will use the frequencies 0,.15,.30,.45...',
		action = FrequencyParser,
		dest = 'frequencies',
		default = [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
	)
	parser.add_argument(
		"-r", "--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into muller_genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint",
		type = float
	)

	parser.add_argument(
		"-l", "--difference-cutoff",
		help = "Minimum p-value to consider a pair of muller_genotypes unrelated. Used when splitting muller_genotypes.",
		action = "store",
		default = 0.25,
		dest = "difference_breakpoint",
		type = float
	)
	parser.add_argument(
		"--genotypes", "--cohorts",
		help = "Indicates that the input table contains muller_genotypes rather than trajectories.",
		action = 'store_true',
		dest = 'is_genotype'
	)
	parser.add_argument(
		"--no-filter",
		help = "Disables genotype filtering.",
		action = 'store_false',
		dest = 'use_filter'
	)
	parser.add_argument(
		"--annotate-all",
		help = "Adds all gene labels to the muller plots, instead of the top three.",
		action = "store_true",
		dest = "annotate_all"
	)
	parser.add_argument(
		"--save-pvalues",
		help = "Saves the p-values to a table and generates a heatmap for the population. Disabling saves a large amount of time for large datasets.",
		action = "store_true",
		dest = "save_pvalue"
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
		default = "matlab",
		dest = "method",
		choices = ['matlab', 'hierarchy']
	)
	parser.add_argument(
		"--metric",
		help = "Selects the distance metric to use. Each metric tends to focus on a specific feature between two series, such as the difference between them or how well they are correlated.",
		action = "store",
		default = "similarity",
		dest = "metric",
		choices = ['similarity', 'binomial', 'pearson', 'minkowski', 'jaccard', 'combined']
	)
	parser.add_argument(
		"-g", "--known-genotypes",
		help = "A file with trajectories known to be in the same genotypes. Each genotype is defined by a comma-delimited line with the labels of the member trajectories.",
		action = "store",
		default = None,
		dest = "known_genotypes"
	)
	parser.add_argument(
		"--genotype-colors",
		help = "An optional map of genotypes to specified colors.",
		action = "store",
		type = Path,
		default = None,
		dest = "genotype_palette_filename"
	)
	return parser


