from loguru import logger
from dataclasses import dataclass

import pandas

from dataio import import_table
from muller.inheritance import order


import pytest


@pytest.fixture
def five_genotypes() -> pandas.DataFrame:
	string = """
	Trajectory	0	1	2	3	5	8	13	21	34
	genotype-C	0	0	0	0.3	0.7	1	1	1	1
	genotype-A	0	0	0	0	0.45	0.5	0.55	0.7	0.85
	genotype-E	0	0	0	0	0	0	0.05	0.55	0.55
	genotype-B	0	0.07	0.1	0.02	0.01	0	0	0	0
	genotype-D	0	0	0	0	0	0	0.07	0	0.01
	"""
	table = import_table(string, index = 'Trajectory')
	return table


@dataclass
class OrderClusterParameters:
	additive_background_double_cutoff: float
	additive_background_single_cutoff: float
	subtractive_background_double_cutoff: float
	subtractive_background_single_cutoff: float
	derivative_detection_cutoff: float
	derivative_check_cutoff: float
	new_background_base_cutoff: float
	new_background_significant_cutoff: float

	@classmethod
	def from_matlab(cls) -> 'OrderClusterParameters':
		return OrderClusterParameters(
			additive_background_double_cutoff = 1.03,
			additive_background_single_cutoff = 1.15,
			subtractive_background_double_cutoff = -0.02,
			subtractive_background_single_cutoff = -0.15,
			derivative_detection_cutoff = 0.02,
			derivative_check_cutoff = 0.01,
			new_background_base_cutoff = 1.0,
			new_background_significant_cutoff = 1.15
		)

	@classmethod
	def from_breakpoints(cls, detection_breakpoint: float, significant_breakpoint: float) -> 'OrderClusterParameters':
		return OrderClusterParameters(
			additive_background_double_cutoff = 1 + detection_breakpoint,
			additive_background_single_cutoff = 1 + significant_breakpoint,
			subtractive_background_double_cutoff = -detection_breakpoint,
			subtractive_background_single_cutoff = -significant_breakpoint,
			derivative_check_cutoff = 0.01,  # Just checking if it's nonnegative
			derivative_detection_cutoff = detection_breakpoint,
			new_background_base_cutoff = 1 + detection_breakpoint,
			new_background_significant_cutoff = 1 + significant_breakpoint
		)


def test_five_genotypes(five_genotypes):
	expected = {
		'genotype-C': 'genotype-0',
		'genotype-A': 'genotype-C',
		'genotype-E': 'genotype-A',
		'genotype-B': 'genotype-0',
		'genotype-D': 'genotype-A'
	}
	options = OrderClusterParameters.from_breakpoints(.03, .15)

	result = order.order_clusters(five_genotypes, options)
	assert result.as_dict() == expected
