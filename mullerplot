#!python
import muller
if __name__ == "__main__":
	args = muller.commandline_parser.create_parser().parse_args()

	workflow = muller.muller_workflow.MullerWorkflow(program_options = args)
	workflow.run(args.filename, args.output_folder)
