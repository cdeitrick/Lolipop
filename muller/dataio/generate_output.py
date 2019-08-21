import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas
from loguru import logger

from muller import dataio, widgets
from muller.clustering.metrics.distance_cache import DistanceCache
from muller.graphics import AnnotatedMullerDiagram, TimeseriesPlot, flowchart, plot_dendrogram, plot_heatmap, palettes

ROOT_GENOTYPE_LABEL = "genotype-0"
FILTERED_GENOTYPE_LABEL = "genotype-filtered"
DEFAULT_GENOTYPE_COLOR = "#333333"


@dataclass
class WorkflowData:
	# Used to organize the output from the workflow.
	version: str
	filename: Path
	program_options: Any
	info: Optional[pandas.DataFrame]
	original_trajectories: Optional[pandas.DataFrame]
	original_genotypes: Optional[pandas.DataFrame]
	trajectories: pandas.DataFrame
	genotypes: pandas.DataFrame
	genotype_members: pandas.Series
	clusters: Any
	p_values: DistanceCache
	filter_cache: List[Tuple[pandas.DataFrame, pandas.DataFrame]]
	linkage_matrix: Any
	genotype_palette_filename: Optional[Path]


class MullerOutputGenerator:
	""" Organizes the workflow for generating the output from the muller workflow."""

	def __init__(self, workflow_data: WorkflowData, output_folder: Path, adjust_populations: bool, render: bool = True):
		self.data = workflow_data
		self.output_folder = output_folder
		self.adjust_populations = adjust_populations
		self.render = render  # Whether to include svg versions of each graph.

		self.timeseries_generator = TimeseriesPlot(render = self.data.program_options['render'])
		self.muller_generator = AnnotatedMullerDiagram(
			outlines = self.data.program_options['draw_outline'],
			render = self.data.program_options['render']
		)
		self.muller_formatter = dataio.GenerateMullerDataFrame()

		self.filenames = dataio.OutputFilenames(output_folder, self.get_base_filename())

		# Read any additional input files
		self.genotype_annotations = self.read_genotype_annotations()
		self.custom_palette = dataio.read_map(self.data.genotype_palette_filename)

		# Now modify the highlighted genotypes
		highlight = self.data.program_options['highlight']

		if highlight:
			highlight_color = self.data.program_options['highlight_color']
			for genotype_identifier in highlight.split(','):
				self.custom_palette[genotype_identifier] = highlight_color

	def run(self) -> None:
		filtered_trajectories = self.save_filtered_trajectories_table()

		edges_table, population_table = self.save_ggmuller_files()
		genotypes: dataio.GenotypeCollection = self.collect_genotype_information(edges_table)

		self.save_base_tables(genotypes.get_trajectory_map())
		self.save_supplementary_files(genotypes)

		logger.info("Generating r script...")
		self.save_r_script(genotypes.get('color_clade'), population_table)

		muller_df = self.muller_formatter.run(edges_table, population_table)
		muller_df.to_csv(self.filenames.table_muller, sep = "\t", index = False)

		logger.info("Generating muller plots...")

		logger.info("Saving linkage files...")
		if self.data.linkage_matrix is not None:
			self.save_linkage_files()

		if self.data.trajectories is not None:
			logger.info("Saving pairwise distances...")
			self.pairwise_distance_information()
		# Save the plots last

		self.save_genotype_plots(genotypes, filtered_trajectories)

		self.save_lineage_plots(genotypes)

		if muller_df is not None:
			self.save_muller_plots(muller_df, genotypes.get('color_clade'), genotypes.get('color_unique'), genotypes.get('annotations'))

	##############################################################################################################################################
	# ---------------------------------------------------------- General Utilites ----------------------------------------------------------------
	##############################################################################################################################################
	def get_base_filename(self) -> str:
		if self.data.program_options['name']:
			base_filename = self.data.program_options['name']
		else:
			base_filename = self.data.filename.stem
		if self.data.program_options['sheetname'] and self.data.program_options['sheetname'] != 'Sheet1':
			base_filename += '.' + self.data.program_options['sheetname']

		return base_filename

	def collect_genotype_information(self, edges_table: pandas.DataFrame) -> dataio.GenotypeCollection:
		""" Consolidates data for each genotype under a single variable. Since many of the attributes of our genotypes end up being
			passed around together, it was a little more convienient to group these attributes in a single dictionary.
		"""

		# Doesn't have to be sorted, but might as well.
		_all_genotype_labels = sorted(set(self.data.original_genotypes.index) | set(self.data.genotypes.index))
		# Generate the color palettes that will be used for all graphics.
		genotype_colors_distinct = palettes.generate_palette(_all_genotype_labels)
		genotype_colors_clade = palettes.generate_palette(edges_table, self.custom_palette, self.genotype_annotations, 'lineage')
		# Hacky method of including the filtered genotype in the palette
		genotype_colors_distinct[FILTERED_GENOTYPE_LABEL] = genotype_colors_clade[FILTERED_GENOTYPE_LABEL] = "#333333"

		genotype_information_collection = dataio.GenotypeCollection()
		for genotype_label in _all_genotype_labels + [FILTERED_GENOTYPE_LABEL, ROOT_GENOTYPE_LABEL]:
			genotype_information = dataio.Genotype(
				label = genotype_label,
				color_custom = None,
				color_clade = genotype_colors_clade[genotype_label],
				color_unique = genotype_colors_distinct[genotype_label],
				annotations = self.genotype_annotations.get(genotype_label, []),
				members = self.data.genotype_members.get(genotype_label, [])
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
	# ---------------------------------------------------------- Generate Output------------------------------------------------------------------
	##############################################################################################################################################
	def save_filtered_trajectories_table(self) -> Optional[pandas.DataFrame]:
		if self.data.original_trajectories is not None:
			filtered_trajectories = self.data.original_trajectories.loc[self.data.genotype_members[FILTERED_GENOTYPE_LABEL]]
			filtered_trajectories['genotype'] = FILTERED_GENOTYPE_LABEL
			merged_table = filtered_trajectories.join(self.data.info, how = 'left')
			merged_table.to_csv(self.filenames.rejected_trajectories, sep = self.filenames.delimiter)
		else:
			merged_table = None
		return merged_table

	def read_genotype_annotations(self) -> Dict[str, List[str]]:
		if self.data.info is not None:
			genotype_annotations = dataio.parse_genotype_annotations(
				self.data.genotype_members,
				self.data.info,
				self.data.program_options['alias_filename']
			)
		else:
			genotype_annotations = {}
		return genotype_annotations

	def save_base_tables(self, parent_genotypes) -> None:
		""" Save the genotype and trajectory tables to the output folder."""
		logger.info("Saving Trajectory and Genotype Tables...")
		self.data.genotypes.to_csv(str(self.filenames.genotype_table), sep = self.filenames.delimiter)

		# Save trajectory tables, if available
		if self.data.trajectories is not None:
			trajectories = dataio.generate_trajectory_table(self.data.trajectories, parent_genotypes, self.data.info)
			trajectories.to_csv(str(self.filenames.trajectory_table), sep = self.filenames.delimiter)

	def save_genotype_plots(self, genotypes: dataio.GenotypeCollection, filtered_trajectories: pandas.DataFrame):
		logger.info("Generating series plots...")
		# Generate the combined trajectory and genotype plots
		# Start by plotting with the clade palette.
		if self.data.trajectories is not None:
			combined_trajectories = pandas.concat([self.data.trajectories, filtered_trajectories])
		else:
			combined_trajectories = None
		self.timeseries_generator.run(
			self.data.genotypes,
			genotypes = genotypes,
			basename = self.filenames.timeseries_plot_genotype_clade,
			palette_type = "color_clade",
			trajectory_timeseries = combined_trajectories
		)
		# The plot by the unique palette
		self.timeseries_generator.run(
			self.data.genotypes,
			genotypes = genotypes,
			basename = self.filenames.timeseries_plot_genotype_unique,
			palette_type = "color_unique",
			trajectory_timeseries = combined_trajectories
		)

		# The inset of the the trajectory plots in the above plots can be a little hard to read.
		# Generate a plot wil only the trajectories. Not sure if its worth including both palette versions.
		if self.data.trajectories is not None:
			trajectory_palette = genotypes.trajectory_palette('color_unique')
			self.timeseries_generator.plot_timeseries(
				timeseries = combined_trajectories,
				palette = trajectory_palette,
				basename = self.filenames.timeseries_plot_trajectory
			)

	def save_ggmuller_files(self) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
		edges_table = self.data.clusters.as_ancestry_table().reset_index()
		edges_table = edges_table[['Parent', 'Identity']]  # Otherwise the r script doesn't work.

		population_table = dataio.generate_ggmuller_population_table(
			self.data.genotypes,
			edges_table,
			self.data.program_options['detection_breakpoint'],
			self.adjust_populations
		)

		population_table.to_csv(str(self.filenames.ggmuller_population), sep = self.filenames.delimiter, index = False)
		edges_table.to_csv(str(self.filenames.ggmuller_edges), sep = self.filenames.delimiter, index = False)
		return edges_table, population_table

	def save_linkage_files(self) -> None:
		num_trajectories = len(self.data.trajectories)
		linkage_table = widgets.format_linkage_matrix(self.data.linkage_matrix, num_trajectories)
		linkage_table.to_csv(str(self.filenames.linkage_matrix_table), sep = self.filenames.delimiter, index = False)
		plot_dendrogram(self.data.linkage_matrix, self.data.p_values, self.filenames.linkage_plot)

	def save_lineage_plots(self, genotypes: dataio.GenotypeCollection) -> None:
		""" Generate the lineage plots."""
		logger.info("Generating Lineage Plots...")
		edges_table = self.data.clusters.priority_table()

		distinct_palette = genotypes.get('color_unique')
		clade_palette = genotypes.get('color_clade')

		genotype_annotations = genotypes.get('annotations')
		flowchart.flowchart(edges_table, clade_palette, annotations = genotype_annotations, filename = self.filenames.lineage_image_clade)
		flowchart.flowchart(edges_table, distinct_palette, annotations = genotype_annotations, filename = self.filenames.lineage_image_distinct)

		if self.render:
			flowchart.flowchart(
				edges_table, clade_palette, annotations = genotype_annotations,
				filename = self.filenames.lineage_image_clade.with_suffix('.svg'))
			flowchart.flowchart(
				edges_table, distinct_palette, annotations = genotype_annotations,
				filename = self.filenames.lineage_image_distinct.with_suffix('.svg'))

	def save_muller_plots(self, muller_df: pandas.DataFrame, clade_palette, distinct_palette, genotype_annotations):
		##############################################################################################################################################
		# ------------------------------------- Generate the muller plot using the table from the r script -------------------------------------------
		##############################################################################################################################################
		logger.info("Generating muller plots...")

		# Draw the muller diagrams
		# Start with the annotated and unannotated clade palettes.
		self.muller_generator.run(muller_df, clade_palette, self.filenames.muller_diagram_clade_annotated, genotype_annotations)
		self.muller_generator.run(muller_df, clade_palette, self.filenames.muller_diagram_clade_unannotated)

		# Draw the distinctive muller diagrams
		self.muller_generator.run(muller_df, distinct_palette, self.filenames.muller_diagram_distinct_annotated, genotype_annotations)
		self.muller_generator.run(muller_df, distinct_palette, self.filenames.muller_diagram_distinct_unannotated)

	def save_r_script(self, palette: Dict[str, str], population_table: pandas.DataFrame) -> None:
		# Generate the rscript and ggmuller DataFrame

		script_content = dataio.generate_r_script(
			trajectory = self.filenames.trajectory_table,
			population = self.filenames.ggmuller_population,
			edges = self.filenames.ggmuller_edges,
			plot_filename = self.filenames.muller_diagram_r_script,
			color_palette = palette,
			genotype_labels = population_table['Identity'].unique().tolist()
		)
		self.filenames.r_script.write_text(script_content)

	def pairwise_distance_information(self) -> None:
		""" Save the pairwise distances for each trajectory."""
		squareform = self.data.p_values.squareform()
		squareform.to_csv(self.filenames.distance_matrix, sep = "\t")

		pvalues_matrix = self.data.p_values.squareform()

		plot_heatmap(pvalues_matrix, self.filenames.distance_heatmap)

	##############################################################################################################################################
	# ---------------------------------------------------- Generate supplementary files ----------------------------------------------------------
	##############################################################################################################################################
	def save_workflow_parameters(self, clade_palette, distinct_palette) -> None:
		parameters = {k: (v if not isinstance(v, Path) else str(v)) for k, v in self.data.program_options.items()}
		parameters['genotypeColorsClade'] = clade_palette
		parameters['genotypeColorsDistinct'] = distinct_palette
		parameters['version'] = self.data.version
		self.filenames.parameters.write_text(json.dumps(parameters, indent = 2))

	def save_supplementary_files(self, genotypes: dataio.GenotypeCollection):
		# Convert the genotypecollection into a regular dictionary so json will play nice.
		_puredict = {k: asdict(v) for k, v in genotypes.items()}
		self.filenames.genotype_information.write_text(json.dumps(_puredict, indent = 2, sort_keys = True))
		self.save_workflow_parameters(genotypes.get('color_clade'), genotypes.get('color_unique'))
		self.data.clusters.to_table().to_csv(self.filenames.lineage_confidence_scores, sep = '\t')
