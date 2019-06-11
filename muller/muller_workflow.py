#!python
"""
	Main script to run the muller workflow.
"""
from pathlib import Path

from loguru import logger

logger.remove()
import sys


try:
	from muller import dataio, clustering, inheritance, commandline_parser
	from muller.muller_output import WorkflowData, generate_output
except ModuleNotFoundError:
	from . import dataio, clustering, inheritance, commandline_parser
	from .muller_output import WorkflowData, generate_output

if commandline_parser.DEBUG:
	logger.add(sys.stderr, level = "INFO")
else:
	logger.add(sys.stderr, level = 'INFO', format="{time:YYYY-MM-DD HH:mm:ss} {level} {message}")

class MullerWorkflow:
	def __init__(self, program_options):
		self.program_options = commandline_parser.parse_workflow_options(program_options)
		logger.info("Program options:")
		for k, v in vars(self.program_options).items():
			logger.info(f"\t{k:<20}{v}")

		breakpoints = self.program_options.frequencies if self.program_options.use_filter else None

		self.genotype_generator = clustering.generate_genotypes.ClusterMutations(
		method = self.program_options.method,
		metric = self.program_options.metric,
		dlimit = self.program_options.detection_breakpoint,
		flimit = self.program_options.fixed_breakpoint,
		sbreakpoint = self.program_options.similarity_breakpoint,
		dbreakpoint = self.program_options.detection_breakpoint,
		breakpoints = breakpoints,
		starting_genotypes = self.program_options.starting_genotypes,
		include_single = self.program_options.include_single
		)

		self.organize_genotypes_workflow = inheritance.SortGenotypeTableWorkflow(
			dlimit = self.program_options.detection_breakpoint,
			slimit = self.program_options.significant_breakpoint,
			flimit = self.program_options.fixed_breakpoint,
			breakpoints = breakpoints
		)

		self.lineage_workflow = inheritance.order.LineageWorkflow(
			dlimit = self.program_options.detection_breakpoint,
			flimit = self.program_options.fixed_breakpoint,
			additive_cutoff = self.program_options.additive_cutoff,
			subtractive_cutoff = self.program_options.subtractive_cutoff,
			derivative_cutoff = self.program_options.derivative_cutoff
		)

	def run(self, filename: Path, output_folder:Path):
		"""
			1. Read input data
			2. Read additional files
			3. calculate the pairwise distances between each trajectory.
			4. Generate Genotypes
			5. Infer lineage
			6. Generate Output.

		"""

		known_ancestry = self.read_additional_files()

		timepoints, mean_genotypes, genotype_members, info = self.generate_genotypes(filename)

		sorted_genotypes = self.organize_genotypes_workflow.run(mean_genotypes)

		genotype_clusters = self.lineage_workflow.run(sorted_genotypes, known_ancestry)

		workflow_data = WorkflowData(
			version = commandline_parser.__VERSION__,
			filename = filename,

			info = info,
			original_trajectories = timepoints,
			original_genotypes = mean_genotypes,
			rejected_trajectories = self.genotype_generator.rejected_trajectories,
			trajectories = timepoints,
			genotypes = sorted_genotypes,
			genotype_members = genotype_members,
			clusters = genotype_clusters,
			program_options = vars(self.program_options),
			p_values = self.genotype_generator.pairwise_distances,
			filter_cache = [],
			linkage_matrix = self.genotype_generator.linkage_table,
			genotype_palette_filename = self.program_options.genotype_palette_filename
		)

		generate_output(
			workflow_data,
			output_folder,
			self.program_options.detection_breakpoint,
			adjust_populations = True
		)


	def generate_genotypes(self, filename:Path):
		if self.program_options.is_genotype:
			mean_genotypes, genotype_info = dataio.parse_genotype_table(filename, self.program_options.sheetname)
			try:
				genotype_members = genotype_info['members']
			except KeyError:
				genotype_members = dict()
			timepoints = info = None

		else:
			timepoints, info = dataio.parse_trajectory_table(filename, self.program_options.sheetname)
			mean_genotypes, genotype_members = self.genotype_generator.run(timepoints)

		return timepoints, mean_genotypes, genotype_members, info

	def read_additional_files(self):
		known_ancestry = dataio.read_map(self.program_options.known_ancestry)
		return known_ancestry




