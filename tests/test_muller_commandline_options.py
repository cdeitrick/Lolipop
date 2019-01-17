import unittest

from muller.commandline_parser import *
from muller.commandline_parser import _parse_frequency_option


class TestMullerProgramOptions(unittest.TestCase):
	def test_parse_commandline_options_set_all_manually(self):
		commandline_parser = create_parser()
		arguments = [
			"--input", str(Path.cwd() / "test_table"),
			"--output", str(Path.cwd() / "output_files"),
			"--fixed", "0.4",
			"--uncertainty", "0.04",
			"--significant", "1.11",
			"-f", "0.3",
			"-r", "0.11",
			"-l", "0.15"
		]

		args = commandline_parser.parse_args(arguments)
		program_options, genotype_options, sort_options, cluster_options = parse_workflow_options(args)
		expected_frequencies = [.9, .6, .3, 0]
		self.assertEqual(Path.cwd() / "test_table", program_options.filename)
		self.assertEqual(Path.cwd() / "output_files", program_options.output_folder)
		self.assertEqual('Sheet1', program_options.sheetname)
		self.assertEqual(.4, program_options.fixed_breakpoint)
		self.assertEqual(.04, program_options.detection_breakpoint)
		self.assertEqual(1.11, program_options.significant_breakpoint)
		self.assertFalse(program_options.mode)
		self.assertListEqual(expected_frequencies, program_options.frequencies)
		self.assertEqual(0.11, program_options.similarity_breakpoint)
		self.assertEqual(0.15, program_options.difference_breakpoint)
		self.assertFalse(program_options.is_genotype)
		self.assertTrue(program_options.use_filter)
		self.assertFalse(program_options.annotate_all)
		self.assertFalse(program_options.annotate_all)

		self.assertEqual(.04, genotype_options.detection_breakpoint)
		self.assertEqual(.4, genotype_options.fixed_breakpoint)
		self.assertIsNone(genotype_options.n_binom)
		self.assertEqual(.11, genotype_options.similarity_breakpoint)
		self.assertEqual(.15, genotype_options.difference_breakpoint)

		self.assertEqual(.04, sort_options.detection_breakpoint)
		self.assertEqual(1.11, sort_options.significant_breakpoint)
		self.assertEqual(.4, sort_options.fixed_breakpoint)
		self.assertListEqual(expected_frequencies, sort_options.frequency_breakpoints)

		self.assertEqual(1 + 1.11, cluster_options.additive_background_single_cutoff)
		self.assertEqual(1.04, cluster_options.additive_background_double_cutoff)
		self.assertEqual(-.04, cluster_options.subtractive_background_double_cutoff)
		self.assertEqual(-1.11, cluster_options.subtractive_background_single_cutoff)
		self.assertEqual(.04, cluster_options.derivative_detection_cutoff)
		self.assertEqual(.01, cluster_options.derivative_check_cutoff)
		self.assertEqual(1.04, cluster_options.new_background_base_cutoff)
		self.assertEqual(1 + 1.11, cluster_options.new_background_significant_cutoff)

	def test_parse_commandline_options_by_breakpoints(self):
		commandline_parser = create_parser()
		detection_cutoff = .012
		arguments = [
			"--input", "test_table",
			"--output", "output_files",
			"--uncertainty", str(detection_cutoff)
		]
		args = commandline_parser.parse_args(arguments)
		program_options, genotype_options, sort_options, cluster_options = parse_workflow_options(args)
		expected_frequencies = [1.0, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0.0]
		self.assertEqual(detection_cutoff, program_options.detection_breakpoint)
		self.assertEqual(1 - detection_cutoff, program_options.fixed_breakpoint)
		self.assertEqual(.15, program_options.significant_breakpoint)
		self.assertListEqual(expected_frequencies, program_options.frequencies)

		self.assertEqual(detection_cutoff, genotype_options.detection_breakpoint)
		self.assertEqual(1 - detection_cutoff, genotype_options.fixed_breakpoint)

		self.assertEqual(detection_cutoff, sort_options.detection_breakpoint)
		self.assertEqual(1 - detection_cutoff, sort_options.fixed_breakpoint)

		self.assertEqual(1 + detection_cutoff, cluster_options.additive_background_double_cutoff)
		self.assertEqual(1.15, cluster_options.additive_background_single_cutoff)
		self.assertEqual(-detection_cutoff, cluster_options.subtractive_background_double_cutoff)
		self.assertEqual(-.15, cluster_options.subtractive_background_single_cutoff)
		self.assertEqual(detection_cutoff, cluster_options.derivative_detection_cutoff)
		self.assertEqual(1 + detection_cutoff, cluster_options.new_background_base_cutoff)
		self.assertEqual(1.15, cluster_options.new_background_significant_cutoff)

	def test_parse_frequency_option(self):
		frequency = "0.3"
		expected_output = [.9, .6, .3, 0]
		output = _parse_frequency_option(frequency)
		self.assertListEqual(expected_output, output)

		frequency = ".2,.4,.7"
		expected_output = [.7, .4, .2, 0]
		self.assertListEqual(expected_output, _parse_frequency_option(frequency))

	def test_parse_lf_options(self):
		commandline_parser = create_parser()
		fixed_cutoff = .5
		difference = .15
		arguments = [
			"--input", "test_table",
			"--output", "output_files",
			"-l", str(difference),
			"--fixed", str(fixed_cutoff)
		]
		args = commandline_parser.parse_args(arguments)
		program_options, genotype_options, sort_options, cluster_options = parse_workflow_options(args)

		self.assertEqual(.03, program_options.detection_breakpoint)
		self.assertEqual(fixed_cutoff, program_options.fixed_breakpoint)
		self.assertEqual(difference, program_options.difference_breakpoint)

if __name__ == "__main__":
	pass