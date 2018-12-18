"""
	Main script to run the muller workflow.
"""
import argparse
from pathlib import Path
from typing import Any, List

from dataclasses import asdict

try:
	from muller.commandline_parser import create_parser, ProgramOptions
	from muller.import_table import import_trajectory_table, import_genotype_table
	from muller_genotypes import calculate_genotypes, sort_genotypes
	from muller import order_clusters
	from muller.muller_output import WorkflowData, generate_output
	from muller import genotype_filters
except ModuleNotFoundError:
	from commandline_parser import create_parser, ProgramOptions
	from import_table import import_trajectory_table, import_genotype_table
	from muller_genotypes import calculate_genotypes, sort_genotypes
	import order_clusters
	import muller_genotypes.sort_genotypes
	import genotype_filters
	from muller_output import WorkflowData, generate_output


def extract_genotypes_from_path(input_filename: Path, sheetname: str):
	mean_genotypes, genotype_info = import_genotype_table(input_filename, sheetname)
	genotype_members = genotype_info['members']
	original_timepoints = timepoints = info = None
	original_genotypes = mean_genotypes
	return original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info


def extract_genotypes_from_trajectories(input_filename: Path, program_options: ProgramOptions, program_options_genotype: Any,
		frequency_breakpoints: List[float]):
	original_timepoints, info = import_trajectory_table(input_filename, program_options.sheetname)
	original_genotypes, genotype_members = calculate_genotypes.workflow(original_timepoints, options = program_options_genotype)

	if program_options.use_filter:
		timepoints, mean_genotypes, genotype_members = genotype_filters.workflow(
			original_timepoints,
			program_options_genotype,
			frequency_breakpoints,
			program_options.use_strict_filter
		)
	else:
		timepoints = original_timepoints.copy()
		mean_genotypes = original_genotypes.copy()

	return original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info

def parse_workflow_options(program_options:argparse.Namespace):
	#program_options = ProgramOptions.from_parser(program_options)
	compatibility_mode = program_options.mode
	if compatibility_mode:
		program_options_genotype = calculate_genotypes.GenotypeOptions.from_matlab()
		program_options_sort = sort_genotypes.SortOptions.from_matlab()
		program_options_clustering = order_clusters.OrderClusterParameters.from_matlab()
	else:
		#freqs = _parse_frequency_option(program_options.frequencies)
		if program_options.fixed_breakpoint is None:
			program_options.fixed_breakpoint = 1 - program_options.detection_breakpoint
		program_options_genotype = calculate_genotypes.GenotypeOptions.from_parser(program_options)
		program_options_clustering = order_clusters.OrderClusterParameters.from_parser(program_options)
		program_options_sort = sort_genotypes.SortOptions(
			detection_breakpoint = program_options_genotype.detection_breakpoint,
			fixed_breakpoint = program_options_genotype.fixed_breakpoint,
			significant_breakpoint = program_options.significant_breakpoint,
			frequency_breakpoints = program_options.frequencies
		)
	return program_options, program_options_genotype, program_options_sort, program_options_clustering

def workflow(input_filename: Path, output_folder: Path, program_options):
	# as long as the sum of the other muller_genotypes that inherit from root is less than 1.
	print("parsing options...")
	program_options, program_options_genotype, program_options_sort, program_options_clustering = parse_workflow_options(program_options)
	from pprint import pprint
	#pprint(asdict(program_options))
	print(program_options)
	print("Importing data...")
	if program_options.is_genotype:
		original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info = extract_genotypes_from_path(input_filename,
			program_options.sheetname)
	else:
		original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info = extract_genotypes_from_trajectories(
			input_filename,
			program_options,
			program_options_genotype,
			program_options_sort.frequency_breakpoints)

	print("sorting muller_genotypes...")
	sorted_genotypes = sort_genotypes.sort_genotypes(mean_genotypes, options = program_options_sort)
	print("nesting muller_genotypes...")
	genotype_clusters = order_clusters.order_clusters(sorted_genotypes, genotype_members, options = program_options_clustering)

	print("Generating output...")
	workflow_data = WorkflowData(
		filename = input_filename,
		info = info,
		original_trajectories = original_timepoints,
		original_genotypes = original_genotypes,
		trajectories = timepoints,
		genotypes = mean_genotypes,
		genotype_members = genotype_members,
		clusters = genotype_clusters,
		genotype_options = program_options_genotype,
		sort_options = program_options_sort,
		cluster_options = program_options_clustering,
		p_values = calculate_genotypes.PAIRWISE_CALCULATIONS,
		filter_cache = []
	)
	generate_output(
		workflow_data,
		output_folder,
		program_options.detection_breakpoint,
		program_options.annotate_all,
		program_options.save_pvalue,
		adjust_populations = True
	)

	return genotype_clusters


if __name__ == "__main__":
	args = create_parser().parse_args()
	#cmd_parser = ProgramOptions.from_parser(args)
	# cmd_parser = ProgramOptions.debug(args)
	workflow(args.filename, args.output_folder, program_options = args)
