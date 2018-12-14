import json
from typing import Any, Dict, Tuple, Optional, List
from pathlib import Path
import pandas
from dataclasses import dataclass

# logging.basicConfig(level = logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')
# logger = logging.getLogger(__name__)
ROOT_GENOTYPE_LABEL = "genotype-0"
FILTERED_GENOTYPE_LABEL = "removed"
OutputType = Tuple[pandas.DataFrame, pandas.DataFrame, str, Dict[str, Any]]

try:
	from muller.muller_genotypes.calculate_genotypes import GenotypeOptions
	from muller_genotypes.sort_genotypes import SortOptions
	from muller.order_clusters import OrderClusterParameters, ClusterType
	from graphics.genotype_plots import plot_genotypes
	from graphics.generate_muller_plot import generate_muller_plot
	from graphics.heatmap import plot_heatmap
	from muller.muller_output.generate_tables import *
	from muller.muller_output.generate_scripts import generate_mermaid_diagram, generate_r_script, excecute_mermaid_script, execute_r_script
	from muller.widgets import generate_genotype_palette, map_trajectories_to_genotype
except ModuleNotFoundError:
	from graphics.genotype_plots import plot_genotypes
	from graphics.generate_muller_plot import generate_muller_plot
	from graphics.heatmap import plot_heatmap
	from muller_output.generate_tables import *
	from muller_output.generate_scripts import generate_mermaid_diagram, generate_r_script, excecute_mermaid_script, execute_r_script
	from widgets import generate_genotype_palette, map_trajectories_to_genotype

	GenotypeOptions = Any
	PairwiseArrayType = Any
	SortOptions = Any
	OrderClusterParameters = Any
	ClusterType = Any


@dataclass
class WorkflowData:
	# Used to organize the output from the workflow.
	filename: Path
	info: Optional[pandas.DataFrame]
	original_trajectories: Optional[pandas.DataFrame]
	original_genotypes: Optional[pandas.DataFrame]
	trajectories: pandas.DataFrame
	genotypes: pandas.DataFrame
	genotype_members: pandas.Series
	clusters: ClusterType
	genotype_options: GenotypeOptions
	sort_options: SortOptions
	cluster_options: OrderClusterParameters
	p_values: PairwiseArrayType
	filter_cache: List[Tuple[pandas.DataFrame, pandas.DataFrame]]


class OutputFilenames:
	# Used to organize the files generated by the workflow.
	def __init__(self, output_folder: Path, name: str, delimiter = '\t'):
		if delimiter == '\t': suffix = 'tsv'
		else: suffix = 'csv'
		subfolder = output_folder / "supplementary-files"
		if not subfolder.exists():
			subfolder.mkdir()
		self.original_trajectory: Path = subfolder / (name + f'.trajectories.original.{suffix}')
		self.original_genotype: Path = subfolder / (name + f'.muller_genotypes.original.{suffix}')
		self.trajectory: Path = output_folder / (name + f'.trajectories.{suffix}')
		self.genotype: Path = output_folder / (name + f'.muller_genotypes.{suffix}')
		self.population: Path = output_folder / (name + f'.ggmuller.populations.{suffix}')
		self.edges: Path = output_folder / (name + f'.ggmuller.edges.{suffix}')
		self.r_script: Path = subfolder / (name + '.r')
		self.muller_table: Path = subfolder / (name + f'.muller.csv')  # This is generated in r.
		self.muller_plot_basic: Path = output_folder / (name + '.muller.basic.png')
		self.muller_plot_annotated: Path = output_folder / (name + '.muller.annotated.png')
		self.mermaid_script: Path = subfolder / (name + '.mermaid.md')
		self.mermaid_render: Path = output_folder / (name + '.mermaid.png')
		self.genotype_plot: Path = output_folder / (name + '.png')
		self.genotype_plot_filtered: Path = output_folder / (name + f".filtered.png")
		self.p_value: Path = subfolder / (name + ".pvalues.tsv")
		self.p_value_matrix: Path = subfolder / (name + f".pvalues.matrix.{suffix}")
		self.p_value_heatmap: Path = subfolder / (name + ".pvalues.heatmap.png")
		self.parameters: Path = output_folder / (name + '.json')


def get_workflow_parameters(workflow_data: WorkflowData, genotype_colors = Dict[str, str]) -> Dict[str, float]:
	parameters = {
		# get_genotype_options
		'detectionCutoff':                        workflow_data.genotype_options.detection_breakpoint,
		'fixedCutoff':                            workflow_data.genotype_options.fixed_breakpoint,
		'similarityCutoff':                       workflow_data.genotype_options.similarity_breakpoint,
		'differenceCutoff':                       workflow_data.genotype_options.difference_breakpoint,
		# sort options
		'significanceCutoff':                     workflow_data.sort_options.significant_breakpoint,
		'frequencyCutoffs':                       workflow_data.sort_options.frequency_breakpoints,
		# cluster options
		'additiveBackgroundDoubleCheckCutoff':    workflow_data.cluster_options.additive_background_double_cutoff,
		'additiveBackgroundSingleCheckCutoff':    workflow_data.cluster_options.additive_background_single_cutoff,
		'subtractiveBackgroundDoubleCheckCutoff': workflow_data.cluster_options.subtractive_background_double_cutoff,
		'subtractiveBackgroundSingleCheckCutoff': workflow_data.cluster_options.subtractive_background_single_cutoff,
		'derivativeDetectionCutoff':              workflow_data.cluster_options.derivative_detection_cutoff,
		'derivativeCheckCutoff':                  workflow_data.cluster_options.derivative_check_cutoff,
		# Palette
		'genotypePalette':                        genotype_colors
	}
	return parameters


def generate_output(workflow_data: WorkflowData, output_folder: Path, detection_cutoff: float, annotate_all: bool, save_pvalues: bool):
	delimiter = '\t'
	parent_genotypes = map_trajectories_to_genotype(workflow_data.genotype_members)

	filtered_trajectories = generate_missing_trajectories_table(workflow_data.trajectories, workflow_data.original_trajectories)
	trajectories = generate_trajectory_table(workflow_data.trajectories, parent_genotypes, workflow_data.info)

	genotype_colors = generate_genotype_palette(workflow_data.original_genotypes.index)
	parameters = get_workflow_parameters(workflow_data, genotype_colors)

	edges_table = generate_ggmuller_edges_table(workflow_data.clusters)
	population_table = generate_ggmuller_population_table(workflow_data.genotypes, edges_table, detection_cutoff)

	filenames = OutputFilenames(output_folder, workflow_data.filename.stem)

	population_table.to_csv(str(filenames.population), sep = delimiter, index = False)
	edges_table.to_csv(str(filenames.edges), sep = delimiter, index = False)

	workflow_data.original_genotypes.to_csv(str(filenames.original_genotype), sep = delimiter)
	workflow_data.genotypes.to_csv(str(filenames.genotype), sep = delimiter)

	if workflow_data.original_trajectories is not None:
		workflow_data.original_trajectories.to_csv(str(filenames.original_trajectory), sep = delimiter)

	if workflow_data.trajectories is not None:
		trajectories.to_csv(str(filenames.trajectory), sep = delimiter)

	muller_df = generate_r_script(
		trajectory = filenames.trajectory,
		population = filenames.population,
		edges = filenames.edges,
		table_filename = filenames.muller_table,
		plot_filename = filenames.muller_plot_basic,
		script_filename = filenames.r_script,
		color_palette = genotype_colors,
		genotype_labels = population_table['Identity'].unique().tolist()
	)
	mermaid_diagram = generate_mermaid_diagram(edges_table, genotype_colors)
	excecute_mermaid_script(filenames.mermaid_script, mermaid_diagram, filenames.mermaid_render)
	filenames.parameters.write_text(json.dumps(parameters, indent = 2))

	plot_genotypes(workflow_data.trajectories, workflow_data.genotypes, filenames.genotype_plot, genotype_colors, parent_genotypes)
	plot_genotypes(filtered_trajectories, workflow_data.genotypes, filenames.genotype_plot_filtered, genotype_colors, parent_genotypes)
	if muller_df is not None:
		generate_muller_plot(muller_df, workflow_data.trajectories, genotype_colors, filenames.muller_plot_annotated, annotate_all)

	if save_pvalues:
		pvalues_table, pvalues_matrix = generate_p_value_table(workflow_data.p_values, parent_genotypes)
		pvalues_table.to_csv(str(filenames.p_value), sep = delimiter, index = False)
		pvalues_matrix.to_csv(str(filenames.p_value_matrix), sep = delimiter, index = False)
		plot_heatmap(pvalues_matrix, filenames.p_value_heatmap)
