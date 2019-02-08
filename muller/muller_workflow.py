"""
	Main script to run the muller workflow.
"""
import logging
from pathlib import Path

logging.basicConfig(filename = "muller_log.txt", level = logging.DEBUG, filemode = 'w', format = '%(module)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__file__)
try:
	from muller.commandline_parser import create_parser, ProgramOptions, parse_workflow_options
	from muller.import_data import import_trajectory_table, import_genotype_table
	from clustering import generate, filters
	from inheritance import order, sort_genotypes
	from muller.muller_output import WorkflowData, generate_output
except ModuleNotFoundError:
	from commandline_parser import create_parser, ProgramOptions, parse_workflow_options
	from import_data import import_trajectory_table, import_genotype_table
	from clustering import generate, filters
	from inheritance import order, sort_genotypes
	from muller_output import WorkflowData, generate_output


def workflow(input_filename: Path, output_folder: Path, program_options):
	# as long as the sum of the other muller_genotypes that inherit from root is less than 1.
	logger.info("parsing options...")
	program_options, program_options_genotype, program_options_sort, program_options_clustering = parse_workflow_options(program_options)
	logger.info("Program options:")
	for k, v in vars(program_options).items():
		logger.info(f"\t{k:<20}{v}")

	logger.info("Importing data...")
	if program_options.is_genotype:
		mean_genotypes, genotype_info = import_genotype_table(input_filename, program_options.sheetname)
		genotype_members = genotype_info['members']
		original_timepoints = timepoints = info = linkage_matrix = None
		original_genotypes = mean_genotypes
	else:
		original_timepoints, info = import_trajectory_table(input_filename, program_options.sheetname)
		original_genotypes, timepoints, mean_genotypes, genotype_members, linkage_matrix = generate.generate_genotypes_with_filter(
			original_timepoints,
			program_options_genotype,
			[program_options.fixed_breakpoint] + program_options.frequencies,
			program_options.use_strict_filter
		)
	logger.info("sorting muller_genotypes...")
	print(mean_genotypes.to_string())
	sorted_genotypes = sort_genotypes.sort_genotypes(mean_genotypes, options = program_options_sort)
	print()
	print(sorted_genotypes.to_string())
	logger.info("nesting muller_genotypes...")
	genotype_clusters = order.order_clusters(sorted_genotypes, options = program_options_clustering)
	logger.info("Generating output...")
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
		p_values = generate.PAIRWISE_CALCULATIONS,
		filter_cache = [],
		linkage_matrix = linkage_matrix
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
	# cmd_parser = ProgramOptions.from_parser(args)
	# cmd_parser = ProgramOptions.debug(args)
	workflow(args.filename, args.output_folder, program_options = args)
