from pathlib import Path


def mullerplot(edges: Path, population: Path, output: Path):
	table_edges = dataio.import_table(edges)
	table_population = dataio.import_table(population)
	from muller.graphics import AnnotatedMullerDiagram
	muller_formatter = dataio.GenerateMullerDataFrame()
	diagram_generator = AnnotatedMullerDiagram()
	muller_df = muller_formatter.run(table_edges, table_population)
	diagram_generator.run(muller_df, output)


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


def main(arguments):
	from muller.muller_workflow import MullerWorkflow
	muller_workflow = MullerWorkflow(arguments)
	muller_workflow.run(arguments.filename, arguments.output_folder)

if __name__ == "__main__":
	import sys
	from pathlib import Path

	sys.path.append(str(Path(__file__).parent.parent))

	from muller import dataio, commandline_parser

	program_parser = commandline_parser.create_parser()
	# Custom method to select `lineage` as the default parser. Used to keep the current api, but will probably be changed later.
	if commandline_parser.DEBUG and False:
		args = [
			"lineage",
			"--sheetname", "BYB1-G07",
			"--input", "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/muller/manuscript/data/nature12344-s2.xlsx",
			"--output", "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/muller/manuscript/data/exampleplots/nature12344",
			"--similarity-cutoff", "0.2"
		]
	else:
		args = sys.argv[1:]

	if args[0] not in {'lineage', 'benchmark', 'muller'}:
		args = ['lineage'] + args
	args = program_parser.parse_args(args)

	if args.name == 'benchmark':
		# The benchmarking utility was activated
		benchmark(args.dataset, args.output)
	elif args.name == 'muller':
		mullerplot(args.edges, args.population, args.output)
	else:
		main(args)

