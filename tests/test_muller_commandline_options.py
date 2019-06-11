from muller.commandline_parser import *
from muller.commandline_parser import _parse_frequency_option


def test_parse_commandline_options_set_all_manually():
	commandline_parser = create_parser()
	arguments = [
		"--input", str(Path.cwd() / "test_table"),
		"--output", str(Path.cwd() / "output_files"),
		"--fixed", "0.4",
		"--detection", "0.04",
		"--significant", "1.11",
		"-f", "0.3",
		"--similarity", "0.11",
		"--difference", "0.15"
	]

	args = commandline_parser.parse_args(arguments)
	program_options = parse_workflow_options(args)
	expected_frequencies = [.9, .6, .3, 0]
	assert program_options.filename == Path.cwd() / "test_table"
	assert program_options.output_folder == Path.cwd() / "output_files"
	assert program_options.sheetname == 0
	assert program_options.fixed_breakpoint == 0.4
	assert program_options.detection_breakpoint == 0.04
	assert program_options.significant_breakpoint == 1.11
	assert program_options.frequencies == expected_frequencies
	assert program_options.similarity_breakpoint == 0.11
	assert program_options.difference_breakpoint == 0.15
	assert not program_options.is_genotype
	assert program_options.use_filter
	assert not program_options.annotate_all


def test_parse_commandline_options_by_breakpoints():
	commandline_parser = create_parser()
	detection_cutoff = .012
	arguments = [
		"--input", "test_table",
		"--output", "output_files",
		"--detection", str(detection_cutoff)
	]
	args = commandline_parser.parse_args(arguments)
	program_options = parse_workflow_options(args)
	expected_frequencies = [1.0, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0.0]
	assert program_options.detection_breakpoint == detection_cutoff
	assert program_options.fixed_breakpoint == 1 - detection_cutoff
	assert program_options.significant_breakpoint == 0.15
	assert program_options.frequencies == expected_frequencies


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
	difference = .15
	arguments = [
		"--input", "test_table",
		"--output", "output_files",
		"--difference-cutoff", str(difference),
		"--fixed", str(fixed_cutoff)
	]
	args = commandline_parser.parse_args(arguments)
	program_options = parse_workflow_options(args)
	assert program_options.detection_breakpoint == 0.03
	assert program_options.fixed_breakpoint == fixed_cutoff
	assert program_options.difference_breakpoint == difference


if __name__ == "__main__":
	pass
