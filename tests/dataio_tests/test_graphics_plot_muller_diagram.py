import pytest
import pandas
from muller import dataio
from muller.graphics import muller_plot
from loguru import logger
@pytest.fixture
def annotation_plot()->muller_plot.AnnotatedMullerDiagram:
	return muller_plot.AnnotatedMullerDiagram(
		outlines = True, render = True
	)

def test_merge_muller_df(annotation_plot):

	index = [0, 8.5, 17, 30.5, 44, 59.5, 75, 90]
	left_series = pandas.DataFrame(
		{
		'Generation': index,
		'Identity': ['genotype-12']*8,
		'Population': [0, 0, 0, 0, 9.4, 4.7, 0, 9.675],
		'Frequency': [0, 0, 0, 0, 0.074, 0.036, 0, 0.044],
		'Group_id': ['genotype-12']*8,
		'Unique_id': [f"genotype-12_{i:.1f}" for i in index]
		}
	)
	right_series = pandas.DataFrame(
		{
		'Generation': index,
		'Identity': ['genotype-12']*8,
		'Population': [0, 0, 0, 0, 9.4, 4.7, 0, 9.675],
		'Frequency': [0, 0, 0, 0, 0.074, 0.036, 0, 0.044],
		'Group_id': ['genotype-12a']*8,
		'Unique_id': [f"genotype-12a_{i:.1f}" for i in index]
		}
	)

	result = annotation_plot._merge_groups(left_series, right_series)
	logger.info(list(result.columns))

	assert list(result['Generation'].values) == index
	assert list(result['Population'].values) == [0, 0, 0, 0, 9.4*2, 4.7*2, 0, 9.675*2]

	assert list(result['Frequency'].values) == [0, 0, 0, 0, 0.074*2, 0.036*2, 0, 0.044*2]
	assert list(result['Identity'].values) == ['genotype-12' for i in range(len(left_series))]
	assert list(result['Group_id'].values) == ['genotype-12' for i in range(len(left_series))]
	assert list(result['Unique_id'].values) == [f'genotype-12_{i:.1f}' for i in index]

def test_merge_muller_table(annotation_plot):
	pass

