from muller.commandline_parser import *
from muller.commandline_parser import _parse_frequency_option


def test_parse_commandline_options_set_all_manually():
	commandline_parser = create_parser()
	arguments = [
		"lineage",
		"--input", str(Path.cwd() / "test_table"),
		"--output", str(Path.cwd() / "output_files"),
		"--fixed", "0.4",
		"--detection", "0.04",
		"--significant", "1.11",
		#"-f", "0.3",
	]

	args = commandline_parser.parse_args(arguments)
	program_options = parse_workflow_options(args)
	expected_frequencies = [.9, .6, .3, 0]
	assert program_options.filename == Path.cwd() / "test_table"
	assert program_options.output_folder == Path.cwd() / "output_files"
	assert program_options.sheetname == 0
	assert program_options.flimit == 0.4
	assert program_options.dlimit == 0.04
	assert program_options.slimit == 1.11
	#assert program_options.frequencies == expected_frequencies
	assert not program_options.is_genotype
	#assert program_options.use_filter
	assert not program_options.annotate_all


def test_parse_commandline_options_by_breakpoints():
	commandline_parser = create_parser()
	detection_cutoff = .012
	arguments = [
		"lineage",
		"--input", "test_table",
		"--output", "output_files",
		"--detection", str(detection_cutoff)
	]
	args = commandline_parser.parse_args(arguments)
	from loguru import logger

	program_options = parse_workflow_options(args)
	for a, b in vars(program_options).items():
		logger.debug(f"{a}\t{b}")
	#expected_frequencies = [1.0, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0.0]
	assert program_options.dlimit == detection_cutoff
	assert program_options.flimit == 1 - detection_cutoff
	assert program_options.slimit == 0.15
	#assert program_options.frequencies == expected_frequencies


def test_parse_frequency_option():
	frequency = "0.3"
	expected_output = [.9, .6, .3, 0]
	output = _parse_frequency_option(frequency)
	assert output == expected_output

	frequency = ".2,.4,.7"
	expected_output = [.7, .4, .2, 0]
	assert _parse_frequency_option(frequency) == expected_output


def test_parse_lf_options():
	commandline_parser = create_parser()
	fixed_cutoff = .5
	arguments = [
		"lineage",
		"--input", "test_table",
		"--output", "output_files",
		"--fixed", str(fixed_cutoff)
	]
	args = commandline_parser.parse_args(arguments)
	program_options = parse_workflow_options(args)
	assert program_options.dlimit == 0.03
	assert program_options.flimit == fixed_cutoff


if __name__ == "__main__":
	pass
