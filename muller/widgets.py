import csv
import itertools
import math
import re
from pathlib import Path
from typing import *

import pandas
from loguru import logger

NumericType = Union[int, float]
IterableValues = Union[List[NumericType], pandas.Series]
NUMERIC_REGEX = re.compile("^.?(?P<number>[\d]+)")

def _coerce_to_series(item:Any)->pandas.Series:
	if not isinstance(item, pandas.Series):
		item = pandas.Series(item)
	return item

def _coerce_to_list(item:Any)->List[Any]:
	if isinstance(item, list):
		# Nothing to change
		result = item
	elif isinstance(item, pandas.Series):
		result = item.tolist()
	else:
		try:
			result = list(item)
		except ValueError:
			message = f"Could not coerce this to a list: '{type(item)}' -> {item}"
			raise ValueError(message)
	return result
def checkdir(path: Union[str, Path]) -> Path:
	path = Path(path)
	if not path.exists():
		path.mkdir()
	return path


def get_numeric_columns(columns: List[Union[str, int, float]]) -> List[Union[str, int, float]]:
	numeric_columns = list()
	for column in columns:
		if isinstance(column, str):
			match = NUMERIC_REGEX.search(column)
			if match:
				col = match.groupdict()['number']
			else:
				continue
		else:
			col = column

		try:
			int(col)
		except (ValueError, TypeError):
			continue
		numeric_columns.append(column)
	return numeric_columns

def get_numeric_table(table:pandas.DataFrame)->pandas.DataFrame:
	""" Returns a version of the input table with only the columns related to timepoints.
		The columns will also be converted to numeric values.
	"""
	table = table[get_numeric_columns(table.column)]
	table.columns = [int(i) for i in table.columns]

	return table

def map_trajectories_to_genotype(genotype_members: pandas.Series) -> Dict[str, str]:
	""" Maps each trajectory to the genotype it belongs to."""
	trajectory_to_genotype = dict()
	for genotype_label, members in genotype_members.items():
		for member in members.split('|'):
			trajectory_to_genotype[member] = genotype_label
	return trajectory_to_genotype


def get_valid_points(left_trajectory: pandas.Series, right_trajectory: pandas.Series, dlimit: float, flimit: Optional[float] = None,
		inner: bool = False) -> Tuple[pandas.Series, pandas.Series]:
	"""
		Filters out timepoints that do not satisfy the detection criteria.
	Parameters
	----------
	left_trajectory: pandas.Series
	right_trajectory: pandas.Series
	dlimit: float
		Removes the timepoint if both points do no exceed this value.
	flimit: float
		If given, all points exceeding this value are treated as "undetected".
	inner: bool; default False
		If `true`, both points must exceed the detection limit for a given timepoint to be considered valid.

	Returns
	-------

	"""

	if flimit is not None:
		# Mask the values so they are considered below the detection limit.
		# Assign a value of -1 so that they will be excluded even if the detection limit is 0.
		# list comprehension is ~6X faster than mask()
		left = pandas.Series([(-1 if s > 0.97 else s) for s in left_trajectory.values], index = left_trajectory.index)
		right = pandas.Series([(-1 if s > 0.97 else s) for s in right_trajectory.values], index = right_trajectory.index)
	else:
		left, right = left_trajectory, right_trajectory

	if inner:
		# list comprehensions are ~10X faster than using the built-in pandas methods.
		at_least_one_detected = [(l > dlimit and r > dlimit) for l, r in zip(left.values, right.values)]
	else:
		at_least_one_detected = [(l > dlimit or r > dlimit) for l, r in zip(left.values, right.values)]
	at_least_one_detected = pandas.Series(at_least_one_detected, index = left.index)

	# Remove indicies where the series value falls below the detection limit. This should include the masked fixed values.
	at_least_one_detected_reduced = at_least_one_detected[at_least_one_detected]
	if at_least_one_detected_reduced.empty:
		# There are no shared timepoints between the series. Assign index_min and index_max to the same number, which will result in an empty dataframe.
		position_index_min = position_index_max = 0
	else:
		# Apparently the min() and max functions now work with strings as well as numbers.
		# Cast the numbers to float so the typeerror is thrown correctly.
		try:
			position_index_min_value = min(at_least_one_detected_reduced.index, key = lambda s: float(s))
			position_index_max_value = max(at_least_one_detected_reduced.index, key = lambda s: float(s))
		except TypeError:
			# The indicies are str and we can't use min() or max(). Assume the indicies are already sorted.
			position_index_min_value = at_least_one_detected_reduced.index[0]
			position_index_max_value = at_least_one_detected_reduced.index[-1]

		position_index_min = left.index.get_loc(position_index_min_value)
		position_index_max = left.index.get_loc(position_index_max_value)
		# Since we want to include the last index, increment position_index_max by one.
		position_index_max += 1
	result_left = left_trajectory[position_index_min:position_index_max]
	result_right = right_trajectory[position_index_min:position_index_max]

	return result_left, result_right


get_detected_points = get_valid_points


def calculate_luminance(color: str) -> float:
	# 0.299 * color.R + 0.587 * color.G + 0.114 * color.B
	red = int(color[1:3], 16)
	green = int(color[3:5], 16)
	blue = int(color[5:], 16)

	lum = (.299 * red) + (.587 * green) + (.114 * blue)
	return lum / 255


def format_inconsistency_matrix(inconsistency_matrix) -> pandas.DataFrame:
	inconsistency_table = pandas.DataFrame(inconsistency_matrix, columns = ['mean', 'std', 'observations', 'statistic'])
	inconsistency_table['observations'] = inconsistency_table['observations'].astype(int)
	return inconsistency_table


# noinspection PyTypeChecker,PyUnresolvedReferences
def fixed_immediately(trajectory: pandas.Series, dlimit: float, flimit: float) -> bool:
	mask1 = (trajectory <= dlimit)
	mask2 = (trajectory >= flimit)
	return (mask1.sum() + mask2.sum()) == len(trajectory)


def fixed(trajectory: pandas.Series, flimit: float) -> bool:
	""" Tests whether the input series fixed at any point.
		If it throws a ValueError due to being an array, check the input data for duplicate labels
	"""
	try:
		return any(i > flimit for i in trajectory.values)
	except ValueError:
		message = "Encountered an error when selecting 'fixed' trajectories. This is usually caused by duplicate trajectory ids in the input table."
		if isinstance(trajectory, pandas.DataFrame):
			message += f"\nRename these trajectories: {list(trajectory.index)}"
		raise ValueError(message)
	except TypeError as exception:
		logger.warning(f"Trajectory: {trajectory.values}")
		logger.warning(f"flimit: {flimit}")

		raise exception

def get_fixed(trajectory:IterableValues, flimit)->pandas.Series:
	""" Returns a pandas.Series object with only the timepoints where the series was fixed."""
	# Usa a dummy value for the `flimit` parameter.
	result = get_intermediate(trajectory, dlimit = flimit, flimit = 2)
	return result
def get_undetected(trajectory:IterableValues, dlimit:float)->pandas.Series:
	return get_intermediate(trajectory, dlimit = -2, flimit = dlimit)
def get_intermediate(trajectory:IterableValues, dlimit:float, flimit:float)->pandas.Series:
	""" Tests whether the input series had timepoints that were neither undetected nor fixed.
		Returns a pandas.Series object with the indecies and values that satisfy this criteria.
	"""
	trajectory = _coerce_to_series(trajectory)
	is_intermediate = trajectory.apply(lambda s: dlimit <= s <= flimit)
	# remove the timepoints which are undetected or fixed
	result = trajectory[is_intermediate]
	return result



def only_fixed(trajectory: pandas.Series, dlimit: float, flimit: float) -> bool:
	""" Tests whether the series immediately fixed and stayed fixed."""
	series = ((i > flimit or i < dlimit) for i in trajectory.values)
	return all(series)

def get_first_fixed_timepoint(elements:pandas.Series, flimit:float)->Any:
	""" Returns the first index that was `fixed` """
	elements = _coerce_to_series(elements)
	result = elements.apply(lambda s: s>flimit)
	# Check if any timepoints were fixed
	if result.sum() == 0:
		result = None
	else:
		result = result.idxmax()

	return result

def get_overlap_regions(left:IterableValues, right:IterableValues, flimit:float)->pandas.Series:
	"""
		Returns the regions between `left` and `right` that have overlapping "fixed" values.
		This will be a pandas.Series object with timepoints as the index and boolean values indicating
		if both series were fixed at that timepoint.
	"""
	# Convert to a pandas.Series object if they aren't already
	left = _coerce_to_series(left)
	right = _coerce_to_series(right)

	fixed_left = left.apply(lambda s: s>= flimit)
	fixed_right = right.apply(lambda s: s>= flimit)
	overlapping_regions = fixed_left & fixed_right
	return overlapping_regions

def overlap(left: Union[List[float], pandas.Series], right: Union[List[float], pandas.Series], dlimit: float) -> int:
	result = [(i > dlimit and j > dlimit) for i, j in zip(left.values, right.values)]

	return sum(result)


# noinspection PyTypeChecker
def find_boundaries_fixed(trajectory: pandas.Series, flimit: float) -> Optional[Tuple[int, int]]:
	fixed_series = trajectory[trajectory > flimit]
	# Assume index is sorted to avoid the situation where the index is str.
	if fixed_series.empty: return None
	else:
		return fixed_series.index[0], fixed_series.index[-1]


def _get_git_log() -> str:
	filename = Path(__file__).parent.parent / ".git" / "logs" / "HEAD"
	try:
		contents = filename.read_text()
	except FileNotFoundError:
		contents = ""
	return contents


def get_commit_hash() -> str:
	commit_hash = "n/a"
	contents = _get_git_log()
	if contents:
		contents = contents.split('\n')
		contents = [i.strip() for i in contents if i.strip()]
		reader = csv.reader(contents, delimiter = '\t')
		for line in reader:
			if line:
				hash_string = line[0]
				try:
					commit_hash = hash_string.split()[1]
				except IndexError:
					continue
		commit_hash = commit_hash[:7]
	else:
		commit_hash = "not available"
	return commit_hash


def get_pair_combinations(elements: Iterable[str]) -> Generator[Tuple[str, str], None, None]:
	# Use a generator to avoid using too much memory for large datasets.
	# The only reason to use a list was so tqdm could automatically figure out how many operations had to be done,
	# but this can be calculated manually by passing the total number of trajectories as an additional parameter.
	return itertools.combinations(elements, 2)


def calculate_number_of_combinations(n: int, r:int = 2) -> int:
	""" Calculates the number of pairwise combinations for an iterable of size `elements` using the
		formula `n!/((n-r)!r!)
		Since we only want combinations of pairs (`r` = 2), this can be simplified to
		`n * (n-1).
	"""
	n_factorial = math.factorial(n)
	r_factorial = math.factorial(r)
	n_minus_r_factorial = math.factorial(n-r)
	#return int(n * (n - r)/r)

	result = n_factorial / (n_minus_r_factorial * r_factorial)
	return int(result)

def validate_table_genotypes(table:pandas.DataFrame):
	""" Throws an error if the genotype table does not conform to the expected format.
		This is mainly for debugging while trying to make sure the genotype table is cosistent at every stage.

		Genotype Table Schema:
		- Columns will be in two forms: integers indicating time points, and a str label for the `members` column.
	"""
	assert 'members' in table.columns
	assert all(isinstance(i, int) for i in table.columns if i != 'members')





def validate_series_edges(edges:pandas.Series):
	""" Validates that the edges table is a pandas.Series object and has correct labels.
		This is ment for debugging while trying to make the edges variable more consistent accross scripts.
	"""

	assert isinstance(edges, pandas.Series)
	assert edges.index.name == 'Identity'
	assert edges.name == 'Parent'