import math

import pandas
import pytest

from muller.inheritance import scoring
from muller.inheritance.genotype_lineage import LineageWorkflow
from .helpers import get_key_pairs, helper_test_score
from ..filenames import real_tables


@pytest.fixture
def scorer() -> scoring.Score:
	# Modify the score weights to make comparing actual vs expected more consistent, since we don't have to worry about ex. 2 vs 1.0 or -1 vs -2.0
	return scoring.Score(0.03, 0.97, 0.05, [1, 1, 1, 1])


@pytest.fixture
def nester() -> LineageWorkflow:
	return LineageWorkflow(0.03, 0.97, 0.05)


@pytest.mark.parametrize("left,right", get_key_pairs(real_tables['nature12344']))
def test_greater_score_nature(scorer, left, right):
	actual_score, expected_score = helper_test_score(scorer, 'greater', left, right, real_tables['nature12344'])

	assert (actual_score == expected_score) or (math.isnan(actual_score) and math.isnan(expected_score))


@pytest.mark.parametrize("left,right", get_key_pairs(real_tables['nature12344']))
def test_score_above_fixed(scorer, left, right):
	actual_score, expected_score = helper_test_score(scorer, 'fixed', left, right, real_tables['nature12344'])
	assert actual_score == expected_score


@pytest.mark.parametrize("left,right", get_key_pairs(real_tables['nature12344']))
def test_score_derivative(scorer, left, right):
	actual_score, expected_score = helper_test_score(scorer, 'derivative', left, right, real_tables['nature12344'])
	assert actual_score == expected_score


@pytest.mark.parametrize("left,right", get_key_pairs(real_tables['nature12344']))
def test_score_jaccard(scorer, left, right):
	actual_score, expected_score = helper_test_score(scorer, 'jaccard', left, right, real_tables['nature12344'])
	assert actual_score == expected_score


def test_lineage():
	nester = LineageWorkflow(0.03, 0.97, 0.05)
	filename = real_tables['nature12344']
	table_genotype = pandas.read_excel(filename, sheet_name = 'genotype').set_index('Genotype')
	expected_lineage = pandas.read_excel(filename, sheet_name = 'edges').set_index('Identity')['Parent']

	actual_lineage = nester.run(table_genotype)

	assert actual_lineage.as_dict() == expected_lineage.to_dict()
