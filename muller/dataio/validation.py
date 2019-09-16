""" Methods to test whether input data is properly formatted. """
import pandas

def test_series_index_match(left:pandas.Series, right:pandas.Series):
	""" Tests whther two series share the same index."""
	assert list(left.index) == list(right.index)

def test_series_index_sorted(series:pandas.Series):
	""" Tests whether the index of a series is properly sorted. """
	assert list(series.index) == sorted(series.index)