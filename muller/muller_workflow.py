"""
	Main script to run the muller workflow.
"""
from pathlib import Path

from loguru import logger

logger.remove()
import sys
logger.add(sys.stderr, level="INFO")
logger.add("muller_log.txt", level = 'DEBUG')
try:
	from muller.commandline_parser import create_parser, parse_workflow_options
	from dataio.trajectories import parse_trajectory_table, parse_genotype_table
	import dataio
	from clustering import generate
	from inheritance import order, sort_genotypes
	from muller.muller_output import WorkflowData, generate_output
except ModuleNotFoundError:
	import dataio
	from commandline_parser import create_parser, parse_workflow_options
	from dataio.trajectories import parse_trajectory_table, parse_genotype_table
	from clustering import generate
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
	if program_options.known_ancestry:
		logger.info(f"Attempting to load ancestry data...")
		known_ancestry = dataio.read_map(program_options.known_ancestry)
	else:
		known_ancestry = dict()

	if program_options.is_genotype:
		mean_genotypes, genotype_info = parse_genotype_table(input_filename, program_options.sheetname)
		try:
			genotype_members = genotype_info['members']
		except KeyError:
			genotype_members = dict()
		timepoints = info = linkage_matrix = None
	else:
		timepoints, info = parse_trajectory_table(input_filename, program_options.sheetname)
		breakpoints = program_options.frequencies if program_options.use_filter else None
		mean_genotypes, genotype_members, linkage_matrix = generate.generate_genotypes(
			timepoints,
			dlimit = program_options.detection_breakpoint,
			flimit = program_options.fixed_breakpoint,
			similarity_breakpoint =  program_options.significant_breakpoint,
			difference_breakpoint =  program_options.difference_breakpoint,
			method = program_options.method,
			metric = program_options.metric,
			breakpoints = breakpoints,
			starting_genotypes = program_options.starting_genotypes
		)

	logger.info("sorting muller_genotypes...")
	sorted_genotypes = sort_genotypes.sort_genotypes(
		mean_genotypes,
		program_options.detection_breakpoint,
		program_options.significant_breakpoint,
		program_options.fixed_breakpoint,
		program_options_sort.frequency_breakpoints
	)
	logger.info("nesting muller_genotypes...")


	genotype_clusters = order.order_clusters(
		sorted_genotypes,
		additive_cutoff = program_options_genotype.detection_breakpoint,
		derivative_cutoff = program_options_clustering.derivative_check_cutoff,
		dlimit = program_options_genotype.detection_breakpoint,
		known_ancestry = known_ancestry
	)
	logger.info("Generating output...")

	workflow_data = WorkflowData(
		filename = input_filename,

		info = info,
		original_trajectories = timepoints,
		original_genotypes = mean_genotypes,
		trajectories = timepoints,
		genotypes = sorted_genotypes,
		genotype_members = genotype_members,
		clusters = genotype_clusters,

		program_options = vars(program_options),
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
		adjust_populations = True
	)
	return genotype_clusters


if __name__ == "__main__":
	args = create_parser().parse_args()

	workflow(args.filename, args.output_folder, program_options = args)
