from pathlib import Path

import pandas
import pytest

from muller import dataio, inheritance, widgets
from muller.inheritance import scoring
from tests.filenames import model_tables
from tests.lineage_tests.helpers import get_key_pairs, helper_test_score
import math
DATA_FOLDER = Path(__file__).parent.parent / "data" / "tables"


############################################################################################################
######################################## Set Up Fixtures ###################################################
############################################################################################################
@pytest.fixture
def scorer() -> scoring.Score:
	# Modify the default weights so we aren't also testing the weights.
	return scoring.Score(0.03, 0.97, 0.05, [1,1,1,1])


@pytest.fixture
def lineage_workflow() -> inheritance.LineageWorkflow:
	runner = inheritance.LineageWorkflow(
		dlimit = 0.03,
		flimit = 0.97,
		pvalue = 0.05
	)
	return runner


@pytest.fixture
def model_periodic_selection() -> pandas.DataFrame:
	filename = model_tables['model.periodicselection']
	table = dataio.import_table(filename, sheet_name = "genotype", index = 'Genotype')
	return table


@pytest.fixture
def model_clonal_interferance() -> pandas.DataFrame:
	filename = model_tables['model.clonalinterferance']
	table = dataio.import_table(filename, sheet_name = 'genotype', index = 'Genotype')
	return table


@pytest.fixture
def model_strong_selection() -> pandas.DataFrame:
	filename = model_tables['model.strongselection']
	table = dataio.import_table(filename, sheet_name = 'genotype', index = "Genotype")
	return table


############################################################################################################
############################################# Helpers ######################################################
############################################################################################################


def helper_for_greater_check(model, left, right):
	left_series = model.loc[left]
	right_series = model.loc[right]

	scores = scoring.Score(0.03, 0.97, 0.05)
	result = scores.calculate_score_greater(left_series, right_series)
	return result


def helper_for_summation_check(model: pandas.DataFrame, left: str, right: str) -> int:
	left_series = model.loc[left]
	right_series = model.loc[right]

	left_series, right_series = widgets.get_valid_points(left_series, right_series, 0.03, inner = False)
	scores = scoring.Score(0.03, 0.97, 0.05)
	result = scores.calculate_score_above_fixed(left_series, right_series)

	return result


def helper_for_area_check(scorer, model: pandas.DataFrame, left: str, right: str):
	left = model.loc[left]
	right = model.loc[right]

	result = scorer.calculate_score_area(left, right)

	return result


############################################################################################################
###################################### Test Derivative Check ###############################################
############################################################################################################
@pytest.mark.skip
@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.strongselection']))
def test_derivative_check_under_strong_selection(scorer,left, right):
	result, expected = helper_test_score(scorer, 'derivative', left, right, model_tables['model.strongselection'])

	assert result == expected


@pytest.mark.skip
@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.periodicselection']))
def test_derivative_check_under_periodic_selection(scorer, left, right):
	actual_score, expected_score = helper_test_score(scorer, 'derivative', left, right, model_tables['model.periodicselection'])
	assert actual_score == expected_score


@pytest.mark.skip
@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.clonalinterferance']))
def test_derivative_check_under_clonal_interferance(scorer, left, right):
	result, expected = helper_test_score(scorer, 'derivative', left, right,model_tables['model.clonalinterferance'] )
	assert result == expected


############################################################################################################
######################################### Test Greater Check ###############################################
############################################################################################################

@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.periodicselection']))
def test_greater_check_under_periodic_selection(scorer, left, right):
	result, expected = helper_test_score(scorer, 'greater', left, right, model_tables['model.periodicselection'])
	assert result == expected


@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.clonalinterferance']))
def test_greater_check_under_clonal_interferance(scorer, left, right):
	result, expected = helper_test_score(scorer, 'greater', left, right, model_tables['model.clonalinterferance'])
	assert result == expected or (math.isnan(result) and math.isnan(expected))


@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.strongselection']))
def test_greater_check_under_strong_selection(scorer, left, right):
	result, expected = helper_test_score(scorer, 'greater', left, right, model_tables['model.strongselection'])
	assert result == expected


############################################################################################################
######################################## Test Above Fixed Check ############################################
############################################################################################################


@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.periodicselection']))
def test_above_fixed_check_under_periodic_selection(scorer, left, right):
	result, expected = helper_test_score(scorer, 'fixed', left, right, model_tables['model.periodicselection'])
	assert result == expected


@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.clonalinterferance']))
def test_above_fixed_check_under_clonal_interferance(scorer, left, right):
	result, expected = helper_test_score(scorer, 'fixed', left, right, model_tables['model.clonalinterferance'])

	assert result == expected


@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.strongselection']))
def test_above_fixed_check_under_strong_selection(scorer, left, right):
	result, expected = helper_test_score(scorer, 'fixed',left, right, model_tables['model.strongselection'])

	assert result == expected


############################################################################################################
############################################ Test Area Check ###############################################
############################################################################################################

@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.periodicselection']))
def test_area_score_periodic_selection(scorer, left, right):
	result, expected = helper_test_score(scorer, 'jaccard', left, right, model_tables['model.periodicselection'])

	assert result == expected


@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.strongselection']))
def test_area_score_strong_selection(scorer, left, right):
	result, expected = helper_test_score(scorer, 'jaccard', left, right, model_tables['model.strongselection'])
	assert result == expected

@pytest.mark.parametrize("left,right", get_key_pairs(model_tables['model.clonalinterferance']))
def test_area_score_clonal_interferance(scorer, left, right):
	result, expected = helper_test_score(scorer, 'jaccard', left, right, model_tables['model.clonalinterferance'])
	assert result == expected

def test_lineage_of_periodic_selection(lineage_workflow, model_periodic_selection):
	expected_lineage = {
		'genotype-red':    'genotype-0',
		'genotype-orange': 'genotype-red',
		'genotype-green':  'genotype-orange',
		'genotype-aqua':   'genotype-green'
	}

	ancestry = lineage_workflow.run(model_periodic_selection).clusters.as_dict()

	assert ancestry == expected_lineage

@pytest.mark.skip
def test_lineage_of_clonal_interferance(lineage_workflow, model_clonal_interferance):
	expected_lineage = {
		'genotype-red':    'genotype-0',
		'genotype-orange': 'genotype-red',
		'genotype-green':  'genotype-red',
		'genotype-orchid': 'genotype-green',
		'genotype-aqua':   'genotype-orange',
		'genotype-sienna': 'genotype-aqua'
	}

	ancestry = lineage_workflow.run(model_clonal_interferance)
	assert ancestry.clusters.as_dict() == expected_lineage


def test_lineage_of_strong_selection(lineage_workflow, model_strong_selection):
	expected_lineage = {
		'genotype-red':    'genotype-0',
		'genotype-orange': 'genotype-red',
		'genotype-green':  'genotype-red',
		'genotype-aqua':   'genotype-green',
		'genotype-orchid': 'genotype-aqua',
		'genotype-sienna': 'genotype-aqua'
	}

	ancestry = lineage_workflow.run(model_strong_selection)
	assert ancestry.clusters.as_dict() == expected_lineage
