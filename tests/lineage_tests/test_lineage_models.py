from pathlib import Path

import pandas
import pytest

from muller import dataio, inheritance, widgets
from muller.inheritance import scoring
from tests.filenames import model_tables

DATA_FOLDER = Path(__file__).parent.parent / "data" / "tables"


############################################################################################################
######################################## Set Up Fixtures ###################################################
############################################################################################################
@pytest.fixture
def scorer() -> scoring.Score:
	return scoring.Score(0.05, 0.03, 0.97)


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
def helper_for_derivative_check(model, left, right) -> bool:
	left_series = model.loc[left]
	right_series = model.loc[right]
	scorer = scoring.Score(0.03, 0.97, 0.05)
	result = scorer.calculate_score_derivative(left_series, right_series)
	return result


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


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		('genotype-red', 'genotype-orange', 2),
		('genotype-orange', 'genotype-green', -2),
		('genotype-orchid', 'genotype-sienna', 0),
		('genotype-orchid', 'genotype-aqua', 2),
		('genotype-red', 'genotype-green', 2),
		('genotype-red', 'genotype-sienna', 0),
		('genotype-orange', 'genotype-orchid', 0)
	]
)
def test_derivative_check_under_strong_selection(model_strong_selection, left, right, expected):
	result = helper_for_derivative_check(model_strong_selection, left, right)

	assert result == expected


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 2),
		("genotype-orange", "genotype-green", 2),
		("genotype-green", "genotype-aqua", 2)
	]
)
def test_derivative_check_under_periodic_selection(model_periodic_selection, left, right, expected):
	result = helper_for_derivative_check(model_periodic_selection, left, right)
	assert result == expected


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 2),
		("genotype-red", "genotype-green", 2),
		("genotype-aqua", "genotype-orchid", 0),
		("genotype-aqua", "genotype-sienna", 2),

		("genotype-orange", "genotype-green", 0),
		# ("genotype-blue", "genotype-orchid", True),
	]
)
def test_derivative_check_under_clonal_interferance(model_clonal_interferance, left, right, expected):
	result = helper_for_derivative_check(model_clonal_interferance, left, right)
	assert result == expected


############################################################################################################
######################################### Test Greater Check ###############################################
############################################################################################################

# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 1),
		("genotype-red", "genotype-green", 1),
		("genotype-red", "genotype-aqua", 1),

		("genotype-orange", "genotype-green", 1),
		("genotype-orange", "genotype-aqua", 1),

		("genotype-green", "genotype-aqua", 1),

		("genotype-orange", "genotype-red", -1),
		("genotype-green", "genotype-orange", -1),
		("genotype-aqua", "genotype-red", -1)

	]
)
def test_greater_check_under_periodic_selection(model_periodic_selection, left, right, expected):
	result = helper_for_greater_check(model_periodic_selection, left, right)

	assert result == expected


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 1),
		("genotype-red", "genotype-green", 1),
		("genotype-red", "genotype-aqua", 1),
		("genotype-red", "genotype-orchid", 1),
		("genotype-red", "genotype-sienna", 1),

		("genotype-orange", "genotype-green", 1),
		("genotype-orange", "genotype-orchid", 1),
		("genotype-orange", "genotype-aqua", 1),
		("genotype-orange", "genotype-sienna", 1),

		("genotype-green", "genotype-orchid", 1),
		("genotype-green", "genotype-aqua", -1),
		("genotype-green", "genotype-sienna", 0),

		("genotype-orchid", "genotype-green", -1),
		("genotype-orchid", "genotype-sienna", -1)
	]
)
def test_greater_check_under_clonal_interferance(model_clonal_interferance, left, right, expected):
	result = helper_for_greater_check(model_clonal_interferance, left, right)
	assert result == expected


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 1),
		("genotype-red", "genotype-green", 1),
		("genotype-red", "genotype-sienna", 1),
		("genotype-red", "genotype-aqua", 1),
		("genotype-red", "genotype-orchid", 1),

		("genotype-orange", "genotype-green", -1),
		("genotype-orange", "genotype-aqua", -1),
		("genotype-orange", "genotype-orchid", -1),
		("genotype-orange", "genotype-sienna", -1),

		("genotype-green", "genotype-orange", -1),
		("genotype-green", "genotype-aqua", 1),
		("genotype-green", "genotype-orchid", 1),
		("genotype-green", "genotype-sienna", 1),

		("genotype-aqua", "genotype-orchid", 1),
		("genotype-aqua", "genotype-sienna", 1),
		("genotype-aqua", "genotype-green", -1),
		("genotype-aqua", "genotype-orange", -1),
		("genotype-aqua", "genotype-red", -1),

		("genotype-orchid", "genotype-aqua", -1),
		("genotype-orchid", "genotype-sienna", -1)
	]
)
def test_greater_check_under_strong_selection(model_strong_selection, left, right, expected):
	result = helper_for_greater_check(model_strong_selection, left, right)
	assert result == expected


############################################################################################################
######################################## Test Above Fixed Check ############################################
############################################################################################################


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 1),
		("genotype-red", "genotype-green", 1),
		("genotype-red", "genotype-aqua", 1),

		("genotype-orange", "genotype-green", 1),
		("genotype-orange", "genotype-aqua", 1),

		("genotype-green", "genotype-aqua", 1),

		("genotype-orange", "genotype-red", 1),
		("genotype-green", "genotype-orange", 1),
		("genotype-aqua", "genotype-red", 1)

	]
)
def test_above_fixed_check_under_periodic_selection(model_periodic_selection, left, right, expected):
	result = helper_for_summation_check(model_periodic_selection, left, right)

	assert result == expected


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 1),
		("genotype-red", "genotype-green", 1),
		("genotype-red", "genotype-orchid", 1),
		("genotype-red", "genotype-aqua", 1),
		("genotype-red", "genotype-sienna", 1),

		("genotype-orange", "genotype-green", 0),
		("genotype-orange", "genotype-orchid", 0),
		("genotype-orange", "genotype-aqua", 1),
		("genotype-orange", "genotype-sienna", 1),

		("genotype-green", "genotype-orchid", 0),
		("genotype-green", "genotype-aqua", 0),
		("genotype-green", "genotype-sienna", 0),

		("genotype-aqua", "genotype-orchid", 0),
		("genotype-aqua", "genotype-sienna", 1)
	]
)
def test_above_fixed_check_under_clonal_interferance(model_clonal_interferance, left, right, expected):
	result = helper_for_summation_check(model_clonal_interferance, left, right)

	assert result == expected


# @pytest.mark.skip
@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 1),
		("genotype-red", "genotype-green", 1),
		("genotype-red", "genotype-sienna", 1),
		("genotype-red", "genotype-aqua", 1),
		("genotype-red", "genotype-orchid", 1),

		("genotype-orange", "genotype-green", 0),
		("genotype-orange", "genotype-aqua", 0),
		("genotype-orange", "genotype-orchid", 0),
		("genotype-orange", "genotype-sienna", 0),

		("genotype-green", "genotype-aqua", 1),
		("genotype-green", "genotype-orchid", 1),
		("genotype-green", "genotype-sienna", 1),

		("genotype-aqua", "genotype-orchid", 0),
		("genotype-aqua", "genotype-sienna", 1),

		("genotype-orchid", "genotype-sienna", 0)
	]
)
def test_above_fixed_check_under_strong_selection(model_strong_selection, left, right, expected):
	result = helper_for_summation_check(model_strong_selection, left, right)

	assert result == expected


############################################################################################################
############################################ Test Area Check ###############################################
############################################################################################################


@pytest.mark.parametrize(
	"left, right, expected",
	[  # nested, unnested
		("genotype-red", "genotype-orange", 2),
		("genotype-red", "genotype-green", 2),
		("genotype-red", "genotype-aqua", 2),

		("genotype-orange", "genotype-green", 2),
		("genotype-orange", "genotype-aqua", 2),

		("genotype-green", "genotype-aqua", 2),

		("genotype-orange", "genotype-red", -2),
		("genotype-green", "genotype-orange", -2),
		("genotype-aqua", "genotype-red", -2)

	]
)
def test_area_score_periodic_selection(scorer, model_periodic_selection, left, right, expected):
	result = helper_for_area_check(scorer, model_periodic_selection, left, right)

	assert result == expected


@pytest.mark.parametrize(
	"left, right, expected",
	[
		("genotype-red", "genotype-orange", 2),
		("genotype-red", "genotype-green", 2),
		("genotype-red", "genotype-sienna", 2),
		("genotype-red", "genotype-aqua", 2),
		("genotype-red", "genotype-orchid", 2),

		("genotype-orange", "genotype-green", -2),
		("genotype-orange", "genotype-aqua", -2),
		("genotype-orange", "genotype-orchid", -2),
		("genotype-orange", "genotype-sienna", -2),

		("genotype-green", "genotype-aqua", 2),
		("genotype-green", "genotype-orchid", 2),
		("genotype-green", "genotype-sienna", 2),

		("genotype-aqua", "genotype-orchid", 2),
		("genotype-aqua", "genotype-sienna", 2),

		("genotype-orchid", "genotype-sienna", -2)
	]
)
def test_area_score_strong_selection(scorer, model_strong_selection, left, right, expected):
	result = helper_for_area_check(scorer, model_strong_selection, left, right)

	assert result == expected


# TODO: Add clonal interference area check


def test_periodic_selection(lineage_workflow, model_periodic_selection):
	expected_lineage = {
		'genotype-red':    'genotype-0',
		'genotype-orange': 'genotype-red',
		'genotype-green':  'genotype-orange',
		'genotype-aqua':   'genotype-green'
	}

	ancestry = lineage_workflow.run(model_periodic_selection).as_dict()

	assert ancestry == expected_lineage


def test_clonal_interferance(lineage_workflow, model_clonal_interferance):
	expected_lineage = {
		'genotype-red':    'genotype-0',
		'genotype-orange': 'genotype-red',
		'genotype-green':  'genotype-red',
		'genotype-orchid': 'genotype-green',
		'genotype-aqua':   'genotype-orange',
		'genotype-sienna': 'genotype-aqua'
	}

	ancestry = lineage_workflow.run(model_clonal_interferance)
	assert ancestry.as_dict() == expected_lineage


def test_strong_selection(lineage_workflow, model_strong_selection):
	expected_lineage = {
		'genotype-red':    'genotype-0',
		'genotype-orange': 'genotype-red',
		'genotype-green':  'genotype-red',
		'genotype-aqua':   'genotype-green',
		'genotype-orchid': 'genotype-aqua',
		'genotype-sienna': 'genotype-aqua'
	}

	ancestry = lineage_workflow.run(model_strong_selection)
	assert ancestry.as_dict() == expected_lineage
