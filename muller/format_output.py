from typing import Any, Dict, Tuple, List, Union, Optional
import pandas
from pathlib import Path
import subprocess
from dataclasses import dataclass

OutputType = Tuple[pandas.DataFrame, pandas.DataFrame, str, Dict[str, Any]]
try:
	import yaml
except ModuleNotFoundError:
	yaml = None

try:
	from muller.get_genotypes import GenotypeOptions
	from muller.sort_genotypes import SortOptions
	from muller.order_clusters import OrderClusterParameters, ClusterType
	from muller.genotype_plots import plot_genotypes
except ModuleNotFoundError:
	# noinspection PyUnresolvedReferences
	from get_genotypes import GenotypeOptions
	from sort_genotypes import SortOptions
	from order_clusters import OrderClusterParameters, ClusterType
	from genotype_plots import plot_genotypes


@dataclass
class WorkflowData:
	filename: Path
	info: Optional[pandas.DataFrame]
	original_trajectories: Optional[pandas.DataFrame]
	original_genotypes: Optional[pandas.DataFrame]
	trajectories: pandas.DataFrame
	genotypes: pandas.DataFrame
	clusters: ClusterType
	genotype_options: GenotypeOptions
	sort_options: SortOptions
	cluster_options: OrderClusterParameters


def get_numeric_columns(columns: List) -> List[Union[int, float]]:
	output_columns = list()
	for col in columns:
		try:
			int(col)
			output_columns.append(col)
		except:
			pass
	return output_columns


def convert_population_to_ggmuller_format(mean_genotypes: pandas.DataFrame) -> pandas.DataFrame:
	table = list()
	# "Generation", "Identity" and "Population"

	numeric_columns = get_numeric_columns(mean_genotypes.columns)
	for index, row in mean_genotypes.iterrows():
		for column, value in row.items():
			# if isinstance(column, str): continue
			if column not in numeric_columns: continue
			line = {
				'Generation': column,
				'Identity':   int(row.name.split('-')[-1]),
				'Population': value * 100
			}
			table.append(line)
	df = pandas.DataFrame(table)
	generations = df.groupby(by = 'Generation')

	root_genotype_generations = list()
	for generation, gdf in generations:
		total = gdf['Population'].sum()
		value = 100.0 - total
		if value < 0:
			value = 0
		row = {
			'Generation': generation,
			'Identity':   0,
			'Population': value
		}
		root_genotype_generations.append(row)

	fdf = pandas.concat([df, pandas.DataFrame(root_genotype_generations)])
	fdf['Identity'] = ["Genotype-{}".format(i) for i in fdf['Identity'].values]
	return fdf.sort_values(by = 'Generation')


def convert_clusters_to_ggmuller_format(mermaid: str) -> pandas.DataFrame:
	lines = mermaid.split('\n')
	table = list()
	for line in lines:
		if '>' not in line: continue
		identity, parent = line.split('-->')
		identity = int(identity.split('-')[-1])
		parent = parent[:-1]
		if parent == 'root':
			parent = 0
		else:
			parent = int(parent.split('-')[-1])
		row = {
			'Parent':   parent,
			'Identity': identity
		}
		table.append(row)
	df = pandas.DataFrame(table)
	df = df[['Parent', 'Identity']]
	df = df.sort_values(by = 'Identity')
	df['Parent'] = ['Genotype-{}'.format(i) for i in df['Parent']]
	df['Identity'] = ['Genotype-{}'.format(i) for i in df['Identity']]
	return df


def generate_mermaid_diagram(clusters: ClusterType) -> str:
	contents = ["graph TD;"]
	for k, cluster in clusters.items():
		if len(cluster.background) > 1:
			cluster_background = cluster.background[::-1][:2]
		elif len(cluster.background) == 1:
			cluster_background = cluster.background + ['root']
		else:
			cluster_background = cluster.background[::-1]

		contents.append("-->".join(cluster_background) + ';')
	return "\n".join(contents)


def generate_r_script(population: Path, edges: Path, output_file: Path, color_palette: List[str]) -> str:
	script = """
	library("ggplot2")
	library("ggmuller")
	
	population <- read.table("{population}", header=TRUE)
	edges <- read.table("{edges}", header=TRUE)
	
	Muller_df <- get_Muller_df(edges, population)
	#Muller_plot(Muller_df, add_legend = TRUE)
	
	palette <- c({palette})
	ggplot(Muller_df, aes_string(x = "Generation", y = "Frequency", group = "Group_id", fill = "Identity", colour = "Identity")) + 
    geom_area() +
    theme(legend.position = "right") +
    guides(linetype = FALSE, color = FALSE) + 
    scale_y_continuous(labels = 25 * (0:4), name = "Percentage", expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) + 
    scale_fill_manual(name = "Identity", values = palette) +
    scale_color_manual(values = palette)
	ggsave("{output}", height = 10, width = 10)
	""".format(
		population = population.absolute(),
		edges = edges.absolute(),
		output = output_file.absolute(),
		palette = ",".join(['"#333333"'] + ['"{}"'.format(i) for i in color_palette])
	)
	script = '\n'.join(i.strip() for i in script.split('\n'))

	return script


def generate_output(workflow_data: WorkflowData, output_folder):
	color_palette = [
		'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
		'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
		'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
		'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080',
		'#ffffff', '#000000'
	]

	population, edges, mermaid, parameters = generate_formatted_output(workflow_data, color_palette)
	save_output(workflow_data, population, edges, mermaid, parameters, output_folder, color_palette)


def get_workflow_parameters(workflow_data: WorkflowData) -> Dict[str, float]:
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
		'derivativeCheckCutoff':                  workflow_data.cluster_options.derivative_check_cutoff
	}
	return parameters


def generate_formatted_output(workflow_data: WorkflowData, color_palette: List[str]) -> OutputType:
	"""
		Generates the final output files for the population.
	Parameters
	----------
	timepoints: pandas.DataFrame
	mean_genotypes: pandas.DataFrame
	clusters: ClusterType
	genotype_options: GenotypeOptions
	sort_options: SortOptions
	cluster_options: OrderClusterParameters

	Returns
	-------
	Tuple[pandas.DataFrame, Dict[str,Any], pandas.DataFrame, Dict[str,float], str, pandas.DataFrame, pandas.DataFrame]
	"""

	workflow_parameters = get_workflow_parameters(workflow_data)
	mermaid_diagram = generate_mermaid_diagram(workflow_data.clusters)

	ggmuller_edge_table = convert_clusters_to_ggmuller_format(mermaid_diagram)
	ggmuller_population_table = convert_population_to_ggmuller_format(workflow_data.genotypes)

	if workflow_data.info is not None:
		genotype_map = {k: v.split('|') for k, v in workflow_data.genotypes['members'].items()}
		trajectory_map = dict()
		for g, v in genotype_map.items():
			for t in v:
				trajectory_map[int(t)] = g
		trajectory_table: pandas.DataFrame = workflow_data.trajectories.copy()

		trajectory_table['genotype'] = [trajectory_map[k] for k in trajectory_table.index]
		trajectory_table = trajectory_table.join(workflow_data.info).sort_values(by = ['genotype'])
		workflow_data.trajectories = trajectory_table

	return ggmuller_population_table, ggmuller_edge_table, mermaid_diagram, workflow_parameters


def save_output(workflow_data: WorkflowData, population_table: pandas.DataFrame, edge_table: pandas.DataFrame,
		mermaid_diagram: str, parameters: Dict, output_folder: Path, color_palette: List[str]):
	name = workflow_data.filename.stem
	delimiter = '\t'
	subfolder = output_folder / "supplementary_files"
	if not subfolder.exists():
		subfolder.mkdir()

	original_trajectory_file = subfolder / (name + '.trajectories.original.tsv')
	trajectory_output_file = output_folder / (name + '.trajectories.tsv')

	original_genotype_file = subfolder / (name + '.genotypes.original.tsv')
	genotype_output_file = output_folder / (name + '.genotypes.tsv')

	population_output_file = output_folder / (name + '.ggmuller.populations.tsv')
	edges_output_file = output_folder / (name + '.ggmuller.edges.tsv')

	r_script_file = subfolder / (name + '.r')
	r_script_graph_file = output_folder / (name + '.muller.png')

	mermaid_diagram_script = subfolder / (name + '.mermaid.md')
	mermaid_diagram_render = output_folder / (name + '.mermaid.png')

	original_genotype_plot_filename = subfolder / (name + '.original.png')
	genotype_plot_filename = output_folder / (name + '.filtered.png')

	population_table.to_csv(str(population_output_file), sep = delimiter, index = False)
	edge_table.to_csv(str(edges_output_file), sep = delimiter, index = False)

	workflow_data.original_genotypes.to_csv(str(original_genotype_file), sep = delimiter)
	workflow_data.genotypes.to_csv(str(genotype_output_file), sep = delimiter)

	if workflow_data.original_trajectories is not None:
		workflow_data.original_trajectories.to_csv(str(original_trajectory_file), sep = delimiter)

	if workflow_data.trajectories is not None:
		workflow_data.trajectories.to_csv(str(trajectory_output_file), sep = delimiter)

	plot_genotypes(workflow_data.original_trajectories, workflow_data.original_genotypes,
		original_genotype_plot_filename, color_palette)
	plot_genotypes(workflow_data.trajectories, workflow_data.genotypes, genotype_plot_filename, color_palette)

	mermaid_diagram_script.write_text(mermaid_diagram)
	try:
		subprocess.call(["mmdc", "-i", mermaid_diagram_script, "-o", mermaid_diagram_render], stdout = subprocess.PIPE,
			stderr = subprocess.PIPE)
	except FileNotFoundError:
		pass

	r_script = generate_r_script(
		population = population_output_file,
		edges = edges_output_file,
		output_file = r_script_graph_file,
		color_palette = color_palette
	)
	r_script_file.write_text(r_script)

	subprocess.call(['Rscript', '--vanilla', '--silent', r_script_file], stdout = subprocess.PIPE,
		stderr = subprocess.PIPE)
	_extra_file = Path.cwd() / "Rplots.pdf"
	if _extra_file.exists():
		_extra_file.unlink()

	if yaml:
		fname = output_folder / (name + '.yaml')
		fname.write_text(yaml.dump(parameters))
	else:
		import json
		fname = output_folder / (name + '.json')
		fname.write_text(json.dumps(parameters, indent = 2))
