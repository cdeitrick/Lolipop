"""
	Main script to run the muller workflow.
"""
import argparse
import itertools
import math
from pathlib import Path
from typing import List, Optional, Union
from dataclasses import dataclass, asdict

try:
	from muller.import_table import import_trajectory_table, import_genotype_table
	from muller_genotypes import calculate_genotypes, sort_genotypes
	from muller import order_clusters
	from muller.muller_output import WorkflowData, generate_output
	from muller import genotype_filters
except ModuleNotFoundError:
	from import_table import import_trajectory_table, import_genotype_table
	from muller_genotypes import calculate_genotypes, sort_genotypes
	import order_clusters
	import muller_genotypes.sort_genotypes
	import genotype_filters
	from muller_output import WorkflowData, generate_output


# For convienience. Helps with autocomplete.
@dataclass
class ProgramOptions:
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
	save_pvalue:bool = True

	def __post_init__(self):
		if self.output_folder and not self.output_folder.exists():
			sup = self.output_folder / "supplementary-files"
			self.output_folder.mkdir()

			if not sup.exists():
				sup.mkdir()

	@classmethod
	def from_parser(cls, parser: Union['ProgramOptions', argparse.Namespace]) -> 'ProgramOptions':
		filename = Path.cwd() / parser.filename
		output_folder = Path.cwd() / parser.output_folder

		fixed = parser.fixed_breakpoint
		if fixed is None:
			fixed = 1 - float(parser.detection_breakpoint)
		else:
			fixed = float(fixed)

		frequencies = _parse_frequency_option(parser.frequencies)

		return cls(
			filename = filename,
			output_folder = output_folder,
			sheetname = parser.sheetname,
			fixed_breakpoint = fixed,
			detection_breakpoint = float(parser.detection_breakpoint),
			significant_breakpoint = float(parser.significant_breakpoint),
			mode = parser.mode,
			frequencies = frequencies,
			similarity_breakpoint = float(parser.similarity_breakpoint),
			difference_breakpoint = float(parser.difference_breakpoint),
			is_genotype = parser.is_genotype,
			use_filter = parser.use_filter,
			annotate_all = parser.annotate_all,
			save_pvalue= parser.save_pvalue
		)

	@classmethod
	def debug(cls, parser) -> 'ProgramOptions':
		parser.filename = "/home/cld100/Documents/github/muller_diagrams/tests/data/B1_Muller.xlsx"
		parser.output_folder = './output'
		parser.use_filter = True
		return cls.from_parser(parser)

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


def parse_workflow_options(program_options):
	compatibility_mode = program_options.mode
	if compatibility_mode:
		program_options_genotype = calculate_genotypes.GenotypeOptions.from_matlab()
		program_options_sort = sort_genotypes.SortOptions.from_matlab()
		program_options_clustering = order_clusters.OrderClusterParameters.from_matlab()
	else:
		freqs = _parse_frequency_option(program_options.frequencies)

		program_options_genotype = calculate_genotypes.GenotypeOptions.from_parser(program_options)
		program_options_clustering = order_clusters.OrderClusterParameters.from_parser(program_options)
		program_options_sort = sort_genotypes.SortOptions(
			detection_breakpoint = program_options_genotype.detection_breakpoint,
			fixed_breakpoint = program_options_genotype.fixed_breakpoint,
			significant_breakpoint = program_options.significant_breakpoint,
			frequency_breakpoints = freqs
		)
	return program_options, program_options_genotype, program_options_sort, program_options_clustering


def workflow(input_filename: Path, output_folder: Path, program_options):
	# as long as the sum of the other muller_genotypes that inherit from root is less than 1.
	print("parsing options...")
	program_options, program_options_genotype, program_options_sort, program_options_clustering = parse_workflow_options(program_options)
	from pprint import pprint
	pprint(asdict(program_options))
	print("Importing data...")
	if program_options.is_genotype:
		mean_genotypes, genotype_info = import_genotype_table(input_filename, program_options.sheetname)
		mean_genotypes['members'] = genotype_info['members']
		original_timepoints = timepoints = info = None
		original_genotypes = mean_genotypes
		filter_cache = list()

	else:
		original_timepoints, info = import_trajectory_table(input_filename, program_options.sheetname)
		original_genotypes = calculate_genotypes.workflow(original_timepoints, options = program_options_genotype)

		if program_options.use_filter:
			timepoints, mean_genotypes, filter_cache = genotype_filters.workflow(
				original_timepoints,
				program_options_genotype,
				program_options_sort.frequency_breakpoints
			)
		else:
			timepoints = original_timepoints.copy()
			mean_genotypes = original_genotypes.copy()
			filter_cache = list()
	print("sorting muller_genotypes...")
	sorted_genotypes = sort_genotypes.workflow(mean_genotypes, options = program_options_sort)
	print("nesting muller_genotypes...")
	genotype_clusters = order_clusters.workflow(sorted_genotypes, options = program_options_clustering)

	print("Generating output...")
	workflow_data = WorkflowData(
		filename = input_filename,
		info = info,
		original_trajectories = original_timepoints,
		original_genotypes = original_genotypes,
		trajectories = timepoints,
		genotypes = mean_genotypes,
		clusters = genotype_clusters,
		genotype_options = program_options_genotype,
		sort_options = program_options_sort,
		cluster_options = program_options_clustering,
		p_values = calculate_genotypes.PAIRWISE_CALCULATIONS,
		filter_cache = filter_cache
	)
	generate_output(
		workflow_data,
		output_folder,
		program_options.detection_breakpoint,
		program_options.annotate_all,
		program_options.save_pvalue
	)

	return genotype_clusters


def create_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename'
	)
	parser.add_argument(
		"--sheetname",
		help = "Indicates the sheet to use if the input table is an excel workbook and the data is not in Sheet1",
		action = 'store',
		dest = 'sheetname',
		default = 'Sheet1'
	)
	parser.add_argument(
		'-o', '--output',
		help = "The folder to save the files to.",
		action = 'store',
		dest = 'output_folder'
	)
	parser.add_argument(
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = 'store',
		dest = 'fixed_breakpoint'
	)
	parser.add_argument(
		"-u", "--uncertainty",
		help = "The uncertainty to apply when performing frequency-based calculations. \
			For example, a frequency at a given timepoint is considered undetected if it falls below 0 + `uncertainty`.",
		action = 'store',
		default = 0.03,
		dest = 'detection_breakpoint'
	)
	parser.add_argument(
		"-s", "--significant",
		help = "The frequency at which to consider a genotype significantly greater than zero.",
		action = 'store',
		default = 0.15,
		dest = "significant_breakpoint"
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
		action = 'store',
		dest = 'frequencies',
		default = 0.10
	)
	parser.add_argument(
		"-r", "--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into muller_genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint"
	)

	parser.add_argument(
		"-l", "--difference-cutoff",
		help = "Minimum p-value to consider a pair of muller_genotypes unrelated. Used when splitting muller_genotypes.",
		action = "store",
		default = 0.10,
		dest = "difference_breakpoint"
	)
	parser.add_argument(
		"--muller_genotypes", "--cohorts",
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

	return parser


if __name__ == "__main__":
	args = create_parser().parse_args()
	cmd_parser = ProgramOptions.from_parser(args)
	#cmd_parser = ProgramOptions.debug(args)
	workflow(cmd_parser.filename, cmd_parser.output_folder, program_options = cmd_parser)
