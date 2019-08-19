if __name__ == "__main__":
	import sys
	from pathlib import Path

	sys.path.append(str(Path(__file__).parent.parent))
	from muller.commandline_parser import create_parser
	from muller.muller_workflow import MullerWorkflow

	args = create_parser().parse_args()
	if args.benchmark:
		from muller.clustering.metrics.calculation_threaded import benchmark, plot_benchmark_results
		_p = Path(__file__)
		filename = _p.parent.parent / "tests/data/tables_input_trajectories/B1_Muller.xlsx"
		benchmark_results = benchmark(filename)
		plot_benchmark_results(benchmark_results)

	else:
		muller_workflow = MullerWorkflow(args)
		muller_workflow.run(args.filename, args.output_folder)
