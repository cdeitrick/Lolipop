from pathlib import Path
from typing import Any, List, Tuple

from commandline_parser import ProgramOptions
from . import filters, generate
try:
	from import_data import import_trajectory_table, import_genotype_table
except ModuleNotFoundError:
	from ..import_data import import_trajectory_table, import_genotype_table


def extract_genotypes(filename: Path, program_options, genotype_options, frequency_breakpoints) -> Tuple:
	if program_options.is_genotype:
		original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info = extract_genotypes_from_path(filename,
			program_options.sheetname)
		linkage_matrix = None
	else:
		original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info, linkage_matrix = extract_genotypes_from_trajectories(
			filename,
			program_options,
			genotype_options,
			frequency_breakpoints)

	return original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info, linkage_matrix


def extract_genotypes_from_path(input_filename: Path, sheetname: str):
	mean_genotypes, genotype_info = import_genotype_table(input_filename, sheetname)
	genotype_members = genotype_info['members']
	original_timepoints = timepoints = info = None
	original_genotypes = mean_genotypes
	return original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info


def extract_genotypes_from_trajectories(input_filename: Path, program_options: ProgramOptions, program_options_genotype: Any,
		frequency_breakpoints: List[float]):
	original_timepoints, info = import_trajectory_table(input_filename, program_options.sheetname)
	original_genotypes, genotype_members, linkage_matrix = generate.generate_genotypes(original_timepoints, options = program_options_genotype)

	if program_options.use_filter:
		timepoints, mean_genotypes, genotype_members, linkage_matrix = filters.filter_genotypes(
			original_timepoints,
			program_options_genotype,
			frequency_breakpoints,
			program_options.use_strict_filter
		)
	else:
		timepoints = original_timepoints.copy()
		mean_genotypes = original_genotypes.copy()

	return original_timepoints, original_genotypes, timepoints, mean_genotypes, genotype_members, info, linkage_matrix
