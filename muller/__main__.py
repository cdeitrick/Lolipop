from pathlib import Path


def mullerplot(edges: Path, population: Path, output: Path):
	table_edges = dataio.import_table(edges)
	table_population = dataio.import_table(population)
	from muller.graphics import AnnotatedMullerDiagram
	muller_formatter = dataio.GenerateMullerDataFrame()
	diagram_generator = AnnotatedMullerDiagram()
	muller_df = muller_formatter.run(table_edges, table_population)
	diagram_generator.run(muller_df, output)


def benchmark(dataset: Path, output: Path):
	from muller.clustering.metrics.distance_calculator import benchmark, plot_benchmark_results
	dataset = dataio.import_table(dataset)
	benchmark_results = benchmark(dataset)
	plot_benchmark_results(benchmark_results, output)


def main(args):
	from muller.muller_workflow import MullerWorkflow
	muller_workflow = MullerWorkflow(args)
	muller_workflow.run(args.filename, args.output_folder)

if __name__ == "__main__":
	import sys
	from pathlib import Path

	sys.path.append(str(Path(__file__).parent.parent))
	from muller.commandline_parser import create_parser

	from muller import dataio

	program_parser = create_parser()
	# Custom method to select `lineage` as the default parser. Used to keep the current api, but will probably be changed later.
	program_parser.set_default_subparser('lineage')  # Only works when using sys.args rather than providing the args directly

	args = program_parser.parse_args()

	if args.name == 'benchmark':
		# The benchmarking utility was activated
		benchmark(args.dataset, args.output)
	elif args.name == 'muller':
		mullerplot(args.edges, args.population, args.output)
	else:
		main(args)

