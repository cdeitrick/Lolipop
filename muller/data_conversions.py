from typing import Any, Dict
import pandas
from pathlib import Path

try:
	from muller.get_genotypes import GenotypeOptions
	from muller.sort_genotypes import SortOptions
	from muller.order_clusters import OrderClusterParameters, ClusterType
except ModuleNotFoundError:
	# noinspection PyUnresolvedReferences
	from get_genotypes import GenotypeOptions
	from sort_genotypes import SortOptions
	from order_clusters import OrderClusterParameters, ClusterType

def convert_population_to_ggmuller_format(mean_genotypes: pandas.DataFrame) -> pandas.DataFrame:
	table = list()
	# "Generation", "Identity" and "Population"
	for index, row in mean_genotypes.iterrows():
		for column, value in row.items():
			if isinstance(column, str):continue

			line = {
				'Generation': column,
				'Identity':   int(row.name.split('-')[-1]),
				'Population': value*100
			}
			table.append(line)
	numeric_columns = [i for i in mean_genotypes.columns if not isinstance(i, str)]
	for column in numeric_columns:
		row = {
			'Generation': column,
			'Identity': 0,
			'Population': 100 if column == min(numeric_columns) else 0
		}
		table.append(row)
	return pandas.DataFrame(sorted(table, key = lambda s: s['Generation']))

def convert_clusters_to_ggmuller_format(mermaid:str)->pandas.DataFrame:
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
			'Parent': parent,
			'Identity': identity
		}
		table.append(row)
	df = pandas.DataFrame(table)
	df = df[['Parent', 'Identity']]
	df = df.sort_values(by = 'Identity')
	return df


def generate_formatted_output(timepoints:pandas.DataFrame, mean_genotypes: pandas.DataFrame, clusters: ClusterType,
		genotype_options: GenotypeOptions, sort_options: SortOptions, cluster_options: OrderClusterParameters) -> Dict[
	str, Any]:

	genotypes = list()

	for label, background in clusters.items():
		genotype_members = label.split('|')
		b = "->".join([i for i in background.background[::-1]])
		genotype = {
			'genotypeLabel':      label,
			'genotypeMembers':    genotype_members,
			'genotypeBackground': b
		}
		genotypes.append(genotype)

	parameters = {
		# get_genotype_options
		'detectionCutoff':                        genotype_options.detection_breakpoint,
		'fixedCutoff':                            genotype_options.fixed_breakpoint,
		'similarityCutoff':                       genotype_options.similarity_breakpoint,
		'differenceCutoff':                       genotype_options.difference_breakpoint,
		# sort options
		'significanceCutoff':                     sort_options.significant_breakpoint,
		'frequencyCutoffs':                       sort_options.frequency_breakpoints,
		# cluster options
		'additiveBackgroundDoubleCheckCutoff':    cluster_options.additive_background_double_cutoff,
		'additiveBackgroundSingleCheckCutoff':    cluster_options.additive_background_single_cutoff,
		'subtractiveBackgroundDoubleCheckCutoff': cluster_options.subtractive_background_double_cutoff,
		'subtractiveBackgroundSingleCheckCutoff': cluster_options.subtractive_background_single_cutoff,
		'derivativeDetectionCutoff':              cluster_options.derivative_detection_cutoff,
		'derivativeCheckCutoff':                  cluster_options.derivative_check_cutoff
	}
	gens = list()
	for label, background in clusters.items():
		genotype_label = background.name
		row = {
			'genotypeLabel':   genotype_label,
			'trajectories':    background.members,
			'timepoints':      list(background.trajectory.index),
			'meanFrequencies': [round(i,4) for i in background.trajectory.tolist()],
			'background':      "->".join([i  for i in background.background[::-1]])
		}
		gens.append(row)


	mermaid_diagram = generate_mermaid_diagram(clusters)
	ggmuller_population_table = convert_population_to_ggmuller_format(mean_genotypes)
	ggmuller_edge_table = convert_clusters_to_ggmuller_format(mermaid_diagram)

	data = {
		'trajectoryTable': timepoints,
		'genotypes':               gens,
		'parameters':              parameters,
		#'genotypeTable':           genotype_data,
		'mermaidDiagram':          mermaid_diagram,
		'ggmullerPopulationTable': ggmuller_population_table,
		'ggmullerEdgeTable': ggmuller_edge_table
	}
	#print(ggmuller_population_table.to_string())
	return data


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

def generate_r_script(population:Path, edges:Path, output_file:Path)->str:

	script = """
	library("ggplot2")
	library("ggmuller")
	
	population <- read.csv("{population}")
	edges <- read.csv("{edges}")
	
	Muller_df <- get_Muller_df(example_edges, example_pop_df)
	Muller_plot(Muller_df)
	
	ggsave("{output}")
	
	""".format(
		population = population.absolute(),
		edges = edges.absolute(),
		output = output_file.absolute()
	)
	script = '\n'.join(i.strip() for i in script.split('\n'))

	return script