#!python
"""
	Main script to run the muller workflow.
"""
from pathlib import Path
from typing import Dict, Tuple

import pandas
from loguru import logger

logger.remove()
import sys

from muller import dataio, clustering, inheritance, commandline_parser
from muller.dataio import generate_output, projectdata

logger.level("COMPLETE", no = 1)
if commandline_parser.DEBUG:
	logger.add(sys.stderr, level = "DEBUG")
else:
	logger.add(sys.stderr, level = 'INFO', format = "{time:YYYY-MM-DD HH:mm:ss} {level} {message}")


class MullerWorkflow:
	def __init__(self, program_options) -> None:
		self.program_options = commandline_parser.parse_workflow_options(program_options)
		logger.info("Program options:")
		for k, v in vars(self.program_options).items():
			logger.info(f"\t{k:<25}{v}")

		breakpoints = self.program_options.frequencies if self.program_options.use_filter else None

		if self.program_options.use_filter:
			self.trajectory_filter = clustering.TrajectoryFilter(
				detection_cutoff = self.program_options.detection_breakpoint,
				fixed_cutoff = self.program_options.fixed_breakpoint,
				filter_consistency = self.program_options.filter_constant,
				filter_single = self.program_options.use_filter_single,
				filter_startfixed = self.program_options.use_filter_startsfixed
			)

		else:
			self.trajectory_filter = None

		self.filter_genotype = clustering.GenotypeFilter(
			detection_cutoff = self.program_options.detection_breakpoint,
			fixed_cutoff = self.program_options.fixed_breakpoint,
			frequencies = self.program_options.frequencies
		)

		self.genotype_generator = clustering.ClusterMutations(
			method = self.program_options.method,
			metric = self.program_options.metric,
			dlimit = self.program_options.detection_breakpoint,
			flimit = self.program_options.fixed_breakpoint,
			pvalue = self.program_options.pvalue,
			dbreakpoint = self.program_options.detection_breakpoint,
			breakpoints = breakpoints,
			starting_genotypes = self.program_options.starting_genotypes,
			trajectory_filter = self.trajectory_filter,
			genotype_filter = self.filter_genotype,
			filename_pairwise = self.program_options.filename_pairwise,
			threads = self.program_options.threads
		)

		self.organize_genotypes_workflow = inheritance.SortGenotypeTableWorkflow(
			dlimit = self.program_options.detection_breakpoint,
			slimit = self.program_options.significant_breakpoint,
			flimit = self.program_options.fixed_breakpoint,
			breakpoints = self.program_options.frequencies
		)

		self.lineage_workflow = inheritance.genotype_lineage.LineageWorkflow(
			dlimit = self.program_options.detection_breakpoint,
			flimit = self.program_options.fixed_breakpoint,
			pvalue = self.program_options.pvalue
		)

	def run(self, filename: Path, output_folder: Path):
		"""
			1. Read input data
			2. Read additional files
			3. calculate the pairwise distances between each trajectory.
			4. Generate Genotypes
			5. Infer lineage
			6. Generate Output.

		"""
		results_basic = projectdata.DataWorkflowBasic(
			version = commandline_parser.__VERSION__,
			filename = filename,
			program_options = vars(self.program_options)
		)

		known_ancestry = self.read_additional_files()

		result_genotype_inference = self.generate_genotypes(filename)

		sorted_genotypes = self.organize_genotypes_workflow.run(result_genotype_inference.table_genotypes)

		results_genotype_lineage, info = self.lineage_workflow.run(sorted_genotypes, known_ancestry)

		# Split up the output data based on its primary function (ex. graphics, data)

		#generate_output.OutputGeneratorFull(results_basic, result_genotype_inference, results_genotype_lineage, workflow_data_graphics).run()

	def generate_genotypes(self, filename: Path) -> Tuple[projectdata.DataGenotypeInference, pandas.DataFrame]:
		if self.program_options.is_genotype:
			mean_genotypes, genotype_info = dataio.parse_genotype_table(filename, self.program_options.sheetname)
			genotype_members = genotype_info.get('members', dict())
			info = None

			output_genotype_inference = projectdata.DataGenotypeInference(
				original_trajectories = None,
				table_genotypes = mean_genotypes,
				genotype_members = genotype_members,
				distance_matrix = None,
				linkage_matrix = None,
				filter_cache = []
			)

		else:
			timepoints, info = dataio.parse_trajectory_table(filename, self.program_options.sheetname)
			output_genotype_inference = self.genotype_generator.run(timepoints)
		return output_genotype_inference, info

	def read_additional_files(self) -> Dict[str, str]:
		known_ancestry = dataio.read_map(self.program_options.known_ancestry)
		return known_ancestry


