if __name__ == "__main__":
	import sys
	from pathlib import Path
	sys.path.append(str(Path(__file__).parent.parent))
	from muller.commandline_parser import create_parser
	from muller.muller_workflow import workflow2
	args = create_parser().parse_args()
	#muller_workflow = MullerWorkflow(args)
	#muller_workflow.run(args.filename, args.output_folder)
	workflow2(args.filename, args.output_folder, program_options = args)

