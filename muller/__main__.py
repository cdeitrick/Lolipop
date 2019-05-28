if __name__ == "__main__":
	import sys
	from pathlib import Path
	sys.path.append(str(Path(__file__).parent.parent))
	from muller.commandline_parser import create_parser
	from muller.muller_workflow import workflow
	args = create_parser().parse_args()

	workflow(args.filename, args.output_folder, program_options = args)