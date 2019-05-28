#!python

from muller import muller_workflow
if __name__ == "__main__":
	args = muller_workflow.create_parser().parse_args()

	muller_workflow.workflow(args.filename, args.output_folder, program_options = args) 
