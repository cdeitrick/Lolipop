import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger

from muller import dataio, graphics
from muller.dataio import projectdata, projectpaths
from muller.graphics import graphicsio

ROOT_GENOTYPE_LABEL = "genotype-0"
FILTERED_GENOTYPE_LABEL = "genotype-filtered"
DEFAULT_GENOTYPE_COLOR = "#333333"


class OutputGeneratorFull:
	""" Organizes the workflow for generating the output from the muller workflow."""

	def __init__(self, data_basic: projectpaths.DataWorkflow, data_genotypes: projectdata.DataGenotypeInference,
			data_linkage: projectdata.DataGenotypeLineage,
			data_graphics: projectdata.DataGraphics):
		adjust_populations: bool = True
		self.render: bool = True  # Whether to include svg versions of each graph.
		self.data_basic = data_basic
		self.data_genotypes = data_genotypes
		self.data_linkage = data_linkage
		self.data_graphics = data_graphics
		self.adjust_populations = adjust_populations
		self.color_filtered_series = "#333333"  # Which color to apply to data series which are missing from the palettes.
		self.paths = dataio.ProjectPaths(data_basic.program_options['output_folder'], self.get_base_filename())

		# Use the graphicio object to generate/save the graphics.
		self.generator_graphics = graphicsio.OutputGeneratorGraphics()

		# Set up the table manipulation workflows.
		self.ggmuller_workflow = dataio.GGMuller(data_basic.program_options['detection_breakpoint'], adjust_populations = self.adjust_populations)
		self.muller_formatter = dataio.GenerateMullerDataFrame()

		# Set up the graphics workflows
		self.timeseries_panel_generator = graphics.TimeseriesPanel(render = self.data_basic.program_options['render'])
		self.timeseries_plot_generator = graphics.TimeseriesPlot(render = self.data_basic.program_options['render'])
		self.muller_generator = graphics.MullerPlot(
			outlines = self.data_basic.program_options['draw_outline'],
			render = self.data_basic.program_options['render']
		)
		self.mullerpanel_generator = graphics.MullerPanel(render = self.data_basic.program_options['render'])

		# Read any additional input files
		self.genotype_annotations = self.read_genotype_annotations()
		self.custom_palette = dataio.read_map(self.data_graphics.genotype_palette_filename)

		# Now modify the highlighted genotypes
		highlight = self.data_basic.program_options['highlight']

		if highlight:
			highlight_color = self.data_basic.program_options['highlight_color']
			for genotype_identifier in highlight.split(','):
				self.custom_palette[genotype_identifier] = highlight_color

	def run(self) -> None:
		# Load the palettes here so they can be included in the supplementary output.
		# Need the resulting filtered trajectories table to include in the graphics.
		self.save_table_genotype()
		self.save_linkage_files()

		filtered_trajectories = self.save_tables_trajectory()
		genotypes, muller_df = self.save_data_ggmuller()
		self.save_supplementary_files(genotypes)
		distance_matrix = self.save_pairwise_distances()

		# Save the plots last
		self.save_graphics(muller_df, genotypes, distance_matrix, filtered_trajectories)

	##############################################################################################################################################
	# ---------------------------------------------------------- General Utilites ----------------------------------------------------------------
	##############################################################################################################################################

	def collect_genotype_information(self, edges_table: pandas.Series) -> dataio.GenotypeCollection:
		""" Consolidates data for each genotype under a single variable. Since many of the attributes of our genotypes end up being
			passed around together, it was a little more convienient to group these attributes in a single dictionary.
		"""

		# Doesn't have to be sorted, but might as well.
		_all_genotype_labels = sorted(set(self.data_genotypes.table_genotypes.index))
		# Generate the color palettes that will be used for all graphics.
		genotype_colors_distinct = graphics.palettes.generate_palette(_all_genotype_labels)
		genotype_colors_clade = graphics.palettes.generate_palette(edges_table.reset_index(), self.custom_palette, self.genotype_annotations,
			'lineage')
		# Hacky method of including the filtered genotype in the palette
		genotype_colors_distinct[FILTERED_GENOTYPE_LABEL] = genotype_colors_clade[FILTERED_GENOTYPE_LABEL] = self.color_filtered_series

		genotype_information_collection = dataio.GenotypeCollection()
		for genotype_label in _all_genotype_labels + [FILTERED_GENOTYPE_LABEL, ROOT_GENOTYPE_LABEL]:
			genotype_information = dataio.Genotype(
				label = genotype_label,
				color_custom = None,
				color_clade = genotype_colors_clade[genotype_label],
				color_unique = genotype_colors_distinct[genotype_label],
				annotations = self.genotype_annotations.get(genotype_label, []),
				members = self.data_genotypes.genotype_members.get(genotype_label, [])
			)
			genotype_information_collection[genotype_label] = genotype_information

		# The custom palette can map either genotype names in the form of "genotype-[\d]+" or in the form of annotations.
		# Read the custom palette
		# Now add the custom color to the corresponding genotype object.
		for genotype_identifier, color in self.custom_palette.items():
			if genotype_identifier.startswith('genotype'):
				selected_genotype = genotype_information_collection[genotype_identifier]
				selected_genotype.color_custom = color
				genotype_information_collection[selected_genotype.label] = selected_genotype
			else:
				selected_genotypes = genotype_information_collection.get_genotype_from_annotation(genotype_identifier)
				if selected_genotypes:
					for genotype_label in selected_genotypes:
						selected_genotype = genotype_information_collection[genotype_label]
						selected_genotype.color_custom = color

		return genotype_information_collection

	##############################################################################################################################################
	# ------------------------------------------------------------ Read input---------------------------------------------------------------------
	##############################################################################################################################################
	def read_genotype_annotations(self) -> Dict[str, List[str]]:
		if self.data_genotypes.info is not None:
			genotype_annotations = dataio.parse_genotype_annotations(
				self.data_genotypes.genotype_members,
				self.data_genotypes.info,
				self.data_basic.program_options['alias_filename']
			)
		else:
			genotype_annotations = {}
		return genotype_annotations

	##############################################################################################################################################
	# ---------------------------------------------------------- Generate Output------------------------------------------------------------------
	##############################################################################################################################################
	def _get_filtered_trajectories(self, original_trajectories: pandas.DataFrame) -> pandas.DataFrame:
		""" Extracts the trajectories that were filtered out by comparing the final trajectories with the input trajectories."""
		filtered_trajectories = original_trajectories.loc[self.data_genotypes.genotype_members[FILTERED_GENOTYPE_LABEL]]
		filtered_trajectories['genotype'] = FILTERED_GENOTYPE_LABEL

		return filtered_trajectories

	@staticmethod
	def _add_genotypes_to_trajectory_table(table_trajectory: pandas.DataFrame, genotype_members: pandas.Series) -> pandas.DataFrame:
		table_trajectory = table_trajectory.copy(deep = True)
		# First, reverse the mapping so that it maps trajectories to parent genotype.
		inverse_members = dict()
		for genotype_label, labels in genotype_members.items():
			for label in labels:
				inverse_members[label] = genotype_label
		# Now add a column to the trajectories table with the relevant genotype.
		table_trajectory['genotype'] = [inverse_members[s] for s in table_trajectory.index]
		return table_trajectory

	def save_tables_trajectory(self) -> Optional[pandas.DataFrame]:
		"""
			Saves the trajectories as two tables: One with all trajectories annotated with genotype, and another with only the filtered trajectories.
		Returns
		-------
		table_trajectories_filtered: pandas.DataFrame
			A table consisting only of trajectories which were filtered out during the filtering step. This should only have the timepoint columns
		"""
		# Unpack the relevant parameters
		# This table should have all trajectories annotated by host genotype.
		table_trajectories_original = self.data_genotypes.original_trajectories
		# Make sure that trajectories were provided.
		if table_trajectories_original is not None:
			table_trajectories_original_annotated = self._add_genotypes_to_trajectory_table(table_trajectories_original,
				self.data_genotypes.genotype_members)
			table_trajectories_original_annotated = table_trajectories_original_annotated.sort_values(by = 'genotype')
			# Save the filtered trajectories in a separate table.
			table_trajectories_filtered = self._get_filtered_trajectories(table_trajectories_original_annotated)
			table_trajectories_filtered_annotated = table_trajectories_filtered.join(self.data_genotypes.info, how = 'left')

			table_trajectories_original.to_csv(self.paths.filename_trajectories_final, sep = self.paths.delimiter)
			table_trajectories_filtered_annotated.to_csv(self.paths.filename_trajectories_filtered, sep = self.paths.delimiter)
		else:
			table_trajectories_filtered = None

		return table_trajectories_filtered

	def save_table_genotype(self) -> None:
		logger.info("Saving Genotype Tables...")
		self.data_genotypes.table_genotypes.to_csv(str(self.paths.filename_table_genotype), sep = self.paths.delimiter)

	def save_data_ggmuller(self) -> Tuple[dataio.GenotypeCollection, pandas.DataFrame]:
		# Get a mapping from child genotype to parent genotype.
		series_edges: pandas.Series = self.data_linkage.clusters.as_ancestry_table()
		table_population = self.ggmuller_workflow.generate_ggmuller_population_table(series_edges, self.data_genotypes.table_genotypes)

		table_population.to_csv(str(self.paths.ggmuller_population), sep = self.paths.delimiter, index = False)
		series_edges.to_csv(str(self.paths.ggmuller_edges), sep = self.paths.delimiter, index = True)  # is pandas.Series, so save the index

		muller_df = self.muller_formatter.run(series_edges, table_population)
		muller_df.to_csv(self.paths.table_muller, sep = "\t", index = False)

		genotypes: dataio.GenotypeCollection = self.collect_genotype_information(series_edges)

		logger.info("Generating r script...")
		self.save_r_script(genotypes.get('color_clade'), table_population)
		return genotypes, muller_df

	def save_linkage_files(self) -> None:
		""" Saves the data files generated when inferring genotype lineage."""
		logger.info("Saving linkage files...")
		if self.data_linkage.score_history:  # Test if it is an empty list
			table_linkage_score_records = pandas.DataFrame(self.data_linkage.score_history)
			table_linkage_score_records.to_csv(self.paths.filename_score_records, sep = self.paths.delimiter, index = False)
		if self.data_linkage.linkage_matrix is not None:
			self.data_linkage.linkage_matrix.to_csv(str(self.paths.linkage_matrix_table), sep = self.paths.delimiter, index = False)

	def save_r_script(self, palette: Dict[str, str], population_table: pandas.DataFrame) -> None:
		# Generate the rscript and ggmuller DataFrame

		script_content = dataio.generate_r_script(
			trajectory = self.paths.trajectory_table,
			population = self.paths.ggmuller_population,
			edges = self.paths.ggmuller_edges,
			plot_filename = self.paths.muller_diagram_r_script,
			color_palette = palette,
			genotype_labels = population_table['Identity'].unique().tolist()
		)
		self.paths.r_script.write_text(script_content)

	def save_pairwise_distances(self) -> pandas.DataFrame:
		""" Save the pairwise distances for each trajectory."""
		squareform = self.data_linkage.p_values.squareform()
		squareform.to_csv(self.paths.distance_matrix, sep = "\t")
		return squareform

	##############################################################################################################################################
	# ---------------------------------------------------- Generate supplementary files ----------------------------------------------------------
	##############################################################################################################################################
	def save_workflow_parameters(self, clade_palette, distinct_palette) -> None:
		parameters = {k: (v if not isinstance(v, Path) else str(v)) for k, v in self.data_basic.program_options.items()}
		parameters['genotypeColorsClade'] = clade_palette
		parameters['genotypeColorsDistinct'] = distinct_palette
		parameters['version'] = self.data_basic.version
		self.paths.parameters.write_text(json.dumps(parameters, indent = 2))

	def save_supplementary_files(self, genotypes: dataio.GenotypeCollection):
		# Convert the genotypecollection into a regular dictionary so json will play nice.
		_puredict = {k: asdict(v) for k, v in genotypes.items()}
		self.paths.genotype_information.write_text(json.dumps(_puredict, indent = 2, sort_keys = True))
		self.save_workflow_parameters(genotypes.get('color_clade'), genotypes.get('color_unique'))
		self.data_linkage.clusters.to_table().to_csv(self.paths.lineage_confidence_scores, sep = '\t')
