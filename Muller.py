#!python
import muller
if __name__ == "__main__":
	args = muller.commandline_parser.create_parser().parse_args()

	muller.muller_workflow.workflow(args.filename, args.output_folder, program_options = args) 
