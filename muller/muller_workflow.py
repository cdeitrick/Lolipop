"""
	Main script to run the muller workflow.
"""
from pathlib import Path

from loguru import logger

logger.remove()
import sys
sys.path.append(str(Path(__file__).parent.parent)) # To deal with import errors.
logger.add(sys.stderr, level = "INFO")
# logger.add("muller_log.txt", level = 'DEBUG')
try:
	from muller.commandline_parser import create_parser, parse_workflow_options
	from muller import dataio, clustering, inheritance
	from muller.muller_output import WorkflowData, generate_output
except ModuleNotFoundError:
	import dataio, clustering, inheritance
	from commandline_parser import create_parser, parse_workflow_options
	from muller_output import WorkflowData, generate_output


def workflow(input_filename: Path, output_folder: Path, program_options):
	# as long as the sum of the other muller_genotypes that inherit from root is less than 1.
	logger.info("parsing options...")
	program_options = parse_workflow_options(program_options)
	logger.info("Program options:")
	for k, v in vars(program_options).items():
		logger.info(f"\t{k:<20}{v}")

	logger.info("Importing data...")
	if program_options.known_ancestry:
		logger.info(f"Attempting to load ancestry data...")
		known_ancestry = dataio.read_map(program_options.known_ancestry)
	else:
		known_ancestry = dict()
	breakpoints = program_options.frequencies if program_options.use_filter else None
	genotype_generator = clustering.generate_genotypes.ClusterMutations(
		method = program_options.method,
		metric = program_options.metric,
		dlimit = program_options.detection_breakpoint,
		flimit = program_options.fixed_breakpoint,
		sbreakpoint = program_options.similarity_breakpoint,
		dbreakpoint = program_options.detection_breakpoint,
		breakpoints = breakpoints,
		starting_genotypes = program_options.starting_genotypes
	)
	if program_options.is_genotype:
		mean_genotypes, genotype_info = dataio.parse_genotype_table(input_filename, program_options.sheetname)
		try:
			genotype_members = genotype_info['members']
		except KeyError:
			genotype_members = dict()
		timepoints = info = None

	else:
		timepoints, info = dataio.parse_trajectory_table(input_filename, program_options.sheetname)
		mean_genotypes, genotype_members = genotype_generator.run(timepoints)


	logger.info("sorting muller_genotypes...")
	sorted_genotypes = inheritance.sort_genotypes(
		mean_genotypes,
		dlimit = program_options.detection_breakpoint,
		slimit = program_options.significant_breakpoint,
		flimit = program_options.fixed_breakpoint,
		breakpoints = program_options.frequencies
	)
	logger.info("nesting muller_genotypes...")

	genotype_clusters = inheritance.order_clusters(
		sorted_genotypes,
		dlimit = program_options.detection_breakpoint,
		flimit = program_options.fixed_breakpoint,
		additive_cutoff = program_options.additive_cutoff,
		subtractive_cutoff = program_options.subtractive_cutoff,
		derivative_cutoff = program_options.derivative_cutoff,
		known_ancestry = known_ancestry
	)
	logger.info("Generating output...")
	# TODO Make 'genotype-0' a variable rather than hard-coding it.
	# TODO Add some options to control how the graphics are generated. Ex. the outlines.
	workflow_data = WorkflowData(
		filename = input_filename,

		info = info,
		original_trajectories = timepoints,
		original_genotypes = mean_genotypes,
		rejected_trajectories = genotype_generator.rejected_trajectories,
		trajectories = timepoints,
		genotypes = sorted_genotypes,
		genotype_members = genotype_members,
		clusters = genotype_clusters,
		program_options = vars(program_options),
		p_values = genotype_generator.pairwise_distances,
		filter_cache = [],
		linkage_matrix = genotype_generator.linkage_table,
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
