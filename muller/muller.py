from pathlib import Path
import argparse
import itertools
import math
from typing import List, Union, Optional
from dataclasses import dataclass

try:
	from muller.order_clusters import ClusterType, OrderClusterParameters
	from muller.get_genotypes import GenotypeOptions
	from muller.sort_genotypes import SortOptions
	from muller.import_table import import_trajectory_table, import_genotype_table

	from muller import get_genotypes
	from muller import order_clusters
	from muller import sort_genotypes
	from muller import format_output
	from muller import genotype_filters
except ModuleNotFoundError:
	# noinspection PyUnresolvedReferences
	from order_clusters import ClusterType, OrderClusterParameters
	# noinspection PyUnresolvedReferences
	from import_table import import_trajectory_table, import_genotype_table
	# noinspection PyUnresolvedReferences
	from get_genotypes import GenotypeOptions
	# noinspection PyUnresolvedReferences
	from sort_genotypes import SortOptions

	# noinspection PyUnresolvedReferences
	import get_genotypes
	# noinspection PyUnresolvedReferences
	import order_clusters
	# noinspection PyUnresolvedReferences
	import sort_genotypes
	# noinspection PyUnresolvedReferences
	import format_output
	# noinspection PyUnresolvedReferences
	import genotype_filters


# For convienience. Helps with autocomplete.


def parse_workflow_options(program_options):
	compatibility_mode = program_options.mode
	if compatibility_mode:
		program_options_genotype = get_genotypes.GenotypeOptions.from_matlab()
		program_options_sort = sort_genotypes.SortOptions.from_matlab()
		program_options_clustering = order_clusters.OrderClusterParameters.from_matlab()
	else:
		freqs = program_options.frequencies
		if isinstance(freqs, str):
			if ',' in freqs:
				freqs = map(float, freqs.split(','))
			else:
				freqs = float(freqs)
		if isinstance(freqs, float):
			freqs = [math.fsum(itertools.repeat(freqs, i)) for i in range(int(1 / freqs) + 1)]
			freqs = [round(i, 2) for i in freqs]

		program_options_genotype = get_genotypes.GenotypeOptions.from_parser(program_options)
		program_options_clustering = order_clusters.OrderClusterParameters.from_parser(program_options)
		program_options_sort = sort_genotypes.SortOptions(
			detection_breakpoint = program_options_genotype.detection_breakpoint,
			fixed_breakpoint = program_options_genotype.fixed_breakpoint,
			significant_breakpoint = program_options.significant_breakpoint,
			frequency_breakpoints = freqs
		)
	return program_options, program_options_genotype, program_options_sort, program_options_clustering


def workflow(input_filename: Path, output_folder: Path, program_options):
	# TODO: The background should be 1-sum(other genotypes) at each timepoint
	# as long as the sum of the other genotypes that inherit from root is less than 1.
	program_options, program_options_genotype, program_options_sort, program_options_clustering = parse_workflow_options(
		program_options)

	if program_options.is_genotype:

		mean_genotypes = import_genotype_table(input_filename, program_options.sheetname)

		original_timepoints = timepoints = info = None
		original_genotypes = mean_genotypes

	else:
		original_timepoints, info = import_trajectory_table(input_filename, program_options.sheetname)

		original_genotypes = get_genotypes.workflow(original_timepoints, options = program_options_genotype)

		if program_options.use_filter:
			timepoints, mean_genotypes = genotype_filters.workflow(input_filename, program_options.detection_breakpoint,
				program_options.fixed_breakpoint)
		else:
			timepoints = original_timepoints.copy()
			mean_genotypes = original_genotypes.copy()

	sorted_genotypes = sort_genotypes.workflow(mean_genotypes)

	genotype_clusters = order_clusters.workflow(sorted_genotypes, options = program_options_clustering)

	workflow_data = format_output.WorkflowData(
		filename = input_filename,
		info = info,
		original_trajectories = original_timepoints,
		original_genotypes = original_genotypes,
		trajectories = timepoints,
		genotypes = mean_genotypes,
		clusters = genotype_clusters,
		genotype_options = program_options_genotype,
		sort_options = program_options_sort,
		cluster_options = program_options_clustering
	)
	format_output.generate_output(workflow_data, output_folder, program_options.fixed_breakpoint)

	return genotype_clusters


@dataclass
class ProgramOptions:
	filename: Path
	output_folder: Path
	sheetname: str = "Sheet1"
	fixed_breakpoint: Optional[float] = None
	detection_breakpoint: float = 0.03
	significant_breakpoint: float = 0.15
	mode: bool = False
	frequencies: Union[float, List[float]] = 0.10
	similarity_breakpoint: float = 0.05
	difference_breakpoint: float = 0.10
	is_genotype: bool = False
	use_filter: bool = True

	def __post_init__(self):
		if self.output_folder and not self.output_folder.exists():
			sup = self.output_folder / "supplementary_files"
			self.output_folder.mkdir()

			if not sup.exists():
				sup.mkdir()

	@classmethod
	def from_parser(cls, parser: Union['ProgramOptions', argparse.Namespace]) -> 'ProgramOptions':
		filename = Path(parser.filename) if parser.filename else None
		output_folder = Path(parser.output_folder) if parser.output_folder else None

		fixed = parser.fixed_breakpoint
		if fixed is None:
			fixed = 1 - parser.detection_breakpoint

		return cls(
			filename = filename,
			output_folder = output_folder,
			sheetname = parser.sheetname,
			fixed_breakpoint = fixed,
			detection_breakpoint = parser.detection_breakpoint,
			significant_breakpoint = parser.significant_breakpoint,
			mode = parser.mode,
			frequencies = parser.frequencies,
			similarity_breakpoint = parser.similarity_breakpoint,
			difference_breakpoint = parser.difference_breakpoint,
			is_genotype = parser.is_genotype,
			use_filter = parser.use_filter
		)

	@classmethod
	def debug(cls, parser)->'ProgramOptions':
		parser.filename = "/home/cld100/Documents/github/muller_diagrams/Data files/B1_muller_try1.xlsx"
		parser.output_folder ='./output'
		parser.use_filter = True
		return cls.from_parser(parser)

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
		help = 'The frequency cutoff to use when sorting the genotypes by first detected frequency. For example, a value of 0.15 will use the frequencies 0,.15,.30,.45...',
		action = 'store',
		dest = 'frequencies',
		default = 0.10
	)
	parser.add_argument(
		"-r", "--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint"
	)

	parser.add_argument(
		"-l", "--difference-cutoff",
		help = "Minimum p-value to consider a pair of genotypes unrelated. Used when splitting genotypes.",
		action = "store",
		default = 0.10,
		dest = "difference_breakpoint"
	)
	parser.add_argument(
		"--genotypes", "--cohorts",
		help = "Indicates that the input table contains genotypes rather than trajectories.",
		action = 'store_true',
		dest = 'is_genotype'
	)
	parser.add_argument(
		"--no-filter",
		help = "Disables genotype filtering.",
		action = 'store_false',
		dest = 'use_filter'
	)
	return parser


if __name__ == "__main__":
	args = create_parser().parse_args()
	cmd_parser = ProgramOptions.from_parser(args)
	DEBUG = False
	if DEBUG:
		# noinspection PyRedeclaration
		cmd_parser = ProgramOptions.debug(args)
	workflow(cmd_parser.filename, cmd_parser.output_folder, program_options = cmd_parser)
