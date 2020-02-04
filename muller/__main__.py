from pathlib import Path


def mullerplot(edges: Path, population: Path, output: Path):
	table_edges = dataio.import_table(edges)
	table_population = dataio.import_table(population)
	from muller.graphics import MullerPlot
	muller_formatter = dataio.GenerateMullerDataFrame()
	diagram_generator = MullerPlot(outlines = True, render = True)
	muller_df = muller_formatter.run(table_edges, table_population)
	diagram_generator.plot(muller_df, output)


def timeseriesplot(filename:Path):
	""" Plots a table of trajectories/genotypes"""
	table = dataio.import_table(filename)
	if 'Trajectory' in table.columns:
		table = table.set_index('Trajectory')
	else:
		table = table.set_index('Genotype')



def benchmark(dataset: Path, output: Path):
	from muller.clustering.metrics.distance_calculator import benchmark, plot_benchmark_results
	dataset = dataio.import_table(dataset)
	benchmark_results = benchmark(dataset)
	plot_benchmark_results(benchmark_results, output)


def main(arguments)->None:
	from muller.workflows.workflow_full import run_workflow
	muller_workflow = run_workflow(arguments)
if __name__ == "__main__":
	import sys
	from pathlib import Path

	sys.path.append(str(Path(__file__).parent.parent))
	from muller import dataio, commandline_parser

	program_parser = commandline_parser.create_parser()
	# Custom method to select `lineage` as the default parser. Used to keep the current api, but will probably be changed later.
	#args = program_parser.parse_args()
	program_arguments = commandline_parser.get_arguments()

	# Need to make sure the argument defaults are applied.

	if program_arguments.name is None:
		program_arguments.name = "lineage"

	if program_arguments.name == 'lineage':
		main(program_arguments)
	else:
		print("Need to use `muller/ lineage --input [input]")

