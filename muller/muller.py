from pathlib import Path
import argparse
from typing import Tuple
import itertools
import math
from typing import Dict
import subprocess

try:
	import yaml
except ModuleNotFoundError:
	yaml = None
try:
	from muller.order_clusters import ClusterType, OrderClusterParameters
	from muller.get_genotypes import GenotypeOptions
	from muller.sort_genotypes import SortOptions
	from muller.import_table import import_trajectory_table, import_genotype_table
	from muller.genotype_plots import plot_genotypes

	from muller import get_genotypes
	from muller import order_clusters
	from muller import sort_genotypes
	from muller import data_conversions
except ModuleNotFoundError:
	from order_clusters import ClusterType, OrderClusterParameters
	from import_table import import_trajectory_table, import_genotype_table
	from get_genotypes import GenotypeOptions
	from sort_genotypes import SortOptions
	from genotype_plots import plot_genotypes
	import get_genotypes
	import order_clusters
	import sort_genotypes
	import data_conversions


def workflow(input_filename: Path, output_folder: Path, program_options):
	# TODO: The background should be 1-sum(other genotypes) at each timepoint
	# as long as the sum of the other genotypes that inherit from root is less than 1.
	compatibility_mode = program_options.mode
	if compatibility_mode:
		program_options_genotype = get_genotypes.GenotypeOptions.from_matlab()
		program_options_sort = sort_genotypes.SortOptions.from_matlab()
		program_options_clustering = order_clusters.OrderClusterParameters.from_matlab()
	else:
		freqs = program_options.frequencies
		if isinstance(freqs, str):
			if ',' in freqs:
				freqs = map(float, freqs.split(','))
			else:
				freqs = float(freqs)
		if isinstance(freqs, float):
			freqs = [math.fsum(itertools.repeat(freqs, i)) for i in range(int(1 / freqs) + 1)]
			freqs = [round(i, 2) for i in freqs]

		program_options_genotype = get_genotypes.GenotypeOptions.from_parser(program_options)
		program_options_clustering = order_clusters.OrderClusterParameters.from_parser(program_options)
		program_options_sort = sort_genotypes.SortOptions(
			detection_breakpoint = program_options_genotype.detection_breakpoint,
			fixed_breakpoint = program_options_genotype.fixed_breakpoint,
			significant_breakpoint = program_options.significant_breakpoint,
			frequency_breakpoints = freqs
		)

	if program_options.is_genotype:
		timepoints = None
		mean_genotypes = import_genotype_table(input_filename)
	else:
		timepoints, info = import_trajectory_table(input_filename)

		mean_genotypes = get_genotypes.workflow(timepoints, options = program_options_genotype)

	sorted_genotypes = sort_genotypes.workflow(mean_genotypes)

	genotype_clusters = order_clusters.workflow(sorted_genotypes, options = program_options_clustering)

	# df = pandas.DataFrame(i.trajectory for i in genotype_clusters.values())
	formatted_output: data_conversions.OutputType = data_conversions.generate_formatted_output(
		timepoints,
		mean_genotypes,
		genotype_clusters,
		program_options_genotype,
		program_options_sort,
		program_options_clustering
	)

	save_output(input_filename, output_folder, formatted_output)
	return genotype_clusters


def save_output(input_file: Path, output_folder: Path, data: data_conversions.OutputType):
	trajectory_table, gens, genotype_table, p, mermaid_diagram, population_table, edge_table = data
	name = input_file.stem
	parameters = {
		'genotypeDescriptions': gens,
		'parameters':           p
	}
	# population_table = data.pop('ggmullerPopulationTable')
	# edge_table = data.pop('ggmullerEdgeTable')
	# genotype_table = data.pop('genotypeTable')
	# trajectory_table = data.pop('trajectoryTable')
	# mermaid_diagram = data.pop('mermaidDiagram')

	population_output_file = output_folder / (name + '.ggmuller_populations.tsv')
	edges_output_file = output_folder / (name + '.ggmuller_edges.tsv')
	genotype_output_file = output_folder / (name + '.genotypes.tsv')
	trajectory_output_file = output_folder / (name + '.trajectories.tsv')
	mermaid_diagram_output = output_folder / (name + '.mermaid')
	r_script_file = output_folder / (name + '.r')
	r_script_graph_file = output_folder / (name + '.muller.png')
	genotype_details_file = output_folder / (name + '.genotypemembers.tsv')
	genotype_plot_filename = output_folder / (name + '.genotypeplot.png')

	plot_genotypes(trajectory_table, genotype_table, genotype_plot_filename)

	with genotype_details_file.open('w') as csv_file:
		for g in gens:
			line = "{}\t".format(g['genotypeLabel']) + "\t".join(g['trajectories'].split('|'))
			csv_file.write(line + '\n')

	population_table.to_csv(str(population_output_file), sep = '\t', index = False)
	edge_table.to_csv(str(edges_output_file), sep = '\t', index = False)
	genotype_table.to_csv(str(genotype_output_file), sep = '\t')
	if trajectory_table is not None:
		trajectory_table.to_csv(str(trajectory_output_file), sep = '\t', index = False)
	mermaid_diagram_output.write_text(mermaid_diagram)

	r_script = data_conversions.generate_r_script(
		population = population_output_file,
		edges = edges_output_file,
		output_file = r_script_graph_file
	)
	r_script_file.write_text(r_script)

	subprocess.call(['Rscript', '--vanilla', '--silent', r_script_file])
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


def create_parser() -> argparse.ArgumentParser:
	parser = argparse.ArgumentParser()
	parser.add_argument(
		'-i', '--input',
		help = "The table of trajectories to cluster.",
		action = 'store',
		dest = 'filename'
	)
	parser.add_argument(
		'-o', '--output',
		help = "The folder to save the files to.",
		action = 'store',
		dest = 'output_folder'
	)
	parser.add_argument(
		'--fixed',
		help = "The minimum frequency at which to consider a mutation fixed.",
		action = 'store',
		dest = 'fixed_breakpoint'
	)
	parser.add_argument(
		"-u", "--uncertainty",
		help = "The uncertainty to apply when performing frequency-based calculations. \
			For example, a frequency at a given timepoint is considered undetected if it falls below 0 + `uncertainty`.",
		action = 'store',
		default = 0.03,
		dest = 'detection_breakpoint'
	)
	parser.add_argument(
		"-s", "--significant",
		help = "The frequency at which to consider a genotype significantly greater than zero.",
		action = 'store',
		default = 0.15,
		dest = "significant_breakpoint"
	)
	parser.add_argument(
		"--matlab",
		help = "Mimics the output of the original matlab script.",
		action = 'store_true',
		dest = "mode"
	)
	parser.add_argument(
		"-f", "--frequencies",
		help = 'The frequency cutoff to use when sorting the genotypes by first detected frequency. For example, a value of 0.15 will use the frequencies 0,.15,.30,.45...',
		action = 'store',
		dest = 'frequencies',
		default = 0.10
	)
	parser.add_argument(
		"-r", "--similarity-cutoff",
		help = "Maximum p-value difference to consider trajectories related. Used when grouping trajectories into genotypes.",
		action = "store",
		default = 0.05,
		dest = "similarity_breakpoint"
	)

	parser.add_argument(
		"-l", "--difference-cutoff",
		help = "Minimum p-value to consider a pair of genotypes unrelated. Used when splitting genotypes.",
		action = "store",
		default = 0.10,
		dest = "difference_breakpoint"
	)
	parser.add_argument(
		"--genotypes", "--cohorts",
		help = "Indicates that the input table contains genotypes rather than trajectories.",
		action = 'store_true',
		dest = 'is_genotype'
	)
	return parser


if __name__ == "__main__":

	cmd_parser = create_parser().parse_args()

	#cmd_parser.is_genotype = True
	#cmd_parser.output_folder = "./B1_muller_try1"

	#cmd_parser.filename = "/home/cld100/Documents/github/muller_diagrams/muller/original/B1_muller_try1.genotypes.tsv"


	_input_filename = Path(cmd_parser.filename)
	if cmd_parser.output_folder:
		_output_folder = Path(cmd_parser.output_folder)
	else:
		_output_folder = _input_filename.parent

	if not _output_folder.exists():
		_output_folder.mkdir()

	workflow(_input_filename, _output_folder, program_options = cmd_parser)
