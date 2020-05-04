from pathlib import Path
from typing import *
import pandas
from loguru import logger
import pytest
def assert_frame_equal(left:pandas.DataFrame, right:pandas.DataFrame):
	""" Tests whether two frames are equal. This is being implemented because the version included
		with pandas checks for things that we don't care about. For example, if dataframe A has an
		index of type `int`, and dataframe B has an index with type `int64`, pandas will fail the test,
		even though these descrepancies are caused because we read in a table rather than using the pipeline to
		generated, so our truthset tables weren't working. More importantly, pandas doesn't have clear debug messages,
		So it will complain that a column is different between the two tables, but won;t report which one.
	"""

	# Make sure the index is equal between the two
	try:
		assert list(left.index) == list(right.index)
	except AssertionError:
		message = f"The index does not match."
		logger.error(message)
		logger.error(f"\t Left index: {list(left.index[:5])}")
		logger.error(f"\t Right index: {list(right.index[:5])}")
		raise AssertionError(message)
	# Make sure the dimensions are the same.
	try:
		assert len(left.index) == len(right.index)
	except AssertionError:
		message = f"The dataframe indecies are not the same size: Left: {len(left.index)}, Right: {len(right.index)}"
		raise AssertionError(message)

	# Make sure both tables have the same columns
	try:
		assert list(left.columns) == list(right.columns)
	except AssertionError:
		message = f"The dataframes have different columns: Left has {len(left.columns)} columns, Right has {len(right.columns)} columns."
		logger.error(message)
		logger.error(f"\tLeft columns: {list(left.columns)}")
		logger.error(f"\tRight columns: {list(right.columns)}")
		raise AssertionError(message)

	# Now iterate over each column and test whether they're equal.

	for column in left.columns:
		l = left[column].tolist()
		r = right[column].tolist()
		is_float = isinstance(l[0], float) and isinstance(r[0], float)

		try:
			# Columns with floating point number need to be approximated.
			if is_float:
				l = [round(i, 3) for i in l]
				r = [round(i, 3) for i in r]
				assert [pytest.approx(i) for i in l] == [pytest.approx(i) for i in r]
			else:
				assert l == r
		except AssertionError:
			message = f"The column '{column}' is not the same between the two tables."
			difference = [i == j for i,j in zip(l, r)]
			logger.error(message)
			logger.error(f"\t Left snapshot ({len(l)} elements of type {l}): {l}")
			logger.error(f"\t Right snapshot ({len(r)} elements of type {r}): {r}")
			logger.error(f"\t Difference: {difference}")
			raise AssertionError(message)

def main():
	pass
if __name__ == "__main__":
	main()