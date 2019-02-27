"""
	Main script to run the muller workflow.
"""
import logging
from pathlib import Path

logging.basicConfig(filename = "muller_log.txt", level = logging.INFO, filemode = 'w', format = '%(module)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__file__)
logger.addHandler(logging.StreamHandler())
try:
	from muller.commandline_parser import create_parser, ProgramOptions, parse_workflow_options
	from dataio.trajectories import parse_trajectory_table, parse_genotype_table
	from clustering import generate, filters
	from inheritance import order, sort_genotypes
	from muller.muller_output import WorkflowData, generate_output
except ModuleNotFoundError:
	from commandline_parser import create_parser, ProgramOptions, parse_workflow_options
	from dataio.trajectories import parse_trajectory_table, parse_genotype_table
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
		mean_genotypes, genotype_info = parse_genotype_table(input_filename, program_options.sheetname)
		try:
			genotype_members = genotype_info['members']
		except KeyError:
			genotype_members = dict()
		original_timepoints = timepoints = info = linkage_matrix = None
		original_genotypes = mean_genotypes
	else:
		original_timepoints, info = parse_trajectory_table(input_filename, program_options.sheetname)
		if program_options.use_filter:
			logger.info("using filter...")
			original_genotypes, timepoints, mean_genotypes, genotype_members, linkage_matrix = generate.generate_genotypes_with_filter(
				original_timepoints,
				program_options_genotype,
				[program_options.fixed_breakpoint] + program_options.frequencies,
				program_options.use_strict_filter
			)
		else:
			logger.info("not using filter...")

			timepoints = original_timepoints
			mean_genotypes, genotype_members, linkage_matrix = generate.generate_genotypes(original_timepoints, program_options_genotype)
			original_genotypes = mean_genotypes
	logger.info("sorting muller_genotypes...")
	sorted_genotypes = sort_genotypes.sort_genotypes(mean_genotypes, options = program_options_sort)
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
		linkage_matrix = linkage_matrix,
		genotype_palette_filename = program_options.genotype_palette_filename
	)
	generate_output(
		workflow_data,
		output_folder,
		program_options.detection_breakpoint,
		program_options.save_pvalue,
		adjust_populations = True
	)

	return genotype_clusters


if __name__ == "__main__":
	args = create_parser().parse_args()
	# cmd_parser = ProgramOptions.from_parser(args)
	# cmd_parser = ProgramOptions.debug(args)
	workflow(args.filename, args.output_folder, program_options = args)
