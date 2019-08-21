if __name__ == "__main__":
	import sys
	from pathlib import Path

	sys.path.append(str(Path(__file__).parent.parent))
	from muller.commandline_parser import create_parser
	from muller.muller_workflow import MullerWorkflow

	program_parser = create_parser()
	# Custom method to select `lineage` as the default parser. Used to keep the current api, but will probably be changed later.
	program_parser.set_default_subparser('lineage')  # Only works when using sys.args rather than providing the args directly
	args = program_parser.parse_args()

	if 'dataset' in dir(args):
		# The benchmarking utility was activated
		from muller.clustering.metrics.distance_calculator import benchmark, plot_benchmark_results

		dataset = args.dataset
		output_filename = args.output
		benchmark_results = benchmark(dataset)
		plot_benchmark_results(benchmark_results, output_filename)
	else:
		muller_workflow = MullerWorkflow(args)
		muller_workflow.run(args.filename, args.output_folder)
