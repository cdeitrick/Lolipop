import csv
import re
from pathlib import Path
from typing import Dict, List, Optional
from loguru import logger
import pandas

NUMERIC_REGEX = re.compile("^.?(?P<number>[\d]+)")


def get_numeric_columns(columns: List[str]) -> List[str]:
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


def map_trajectories_to_genotype(genotype_members: pandas.Series) -> Dict[str, str]:
	""" Maps each trajectory to the genotype it belongs to."""
	trajectory_to_genotype = dict()
	for genotype_label, members in genotype_members.items():
		for member in members.split('|'):
			trajectory_to_genotype[member] = genotype_label
	return trajectory_to_genotype


def get_valid_points(left: pandas.Series, right: pandas.Series, dlimit: float, flimit: Optional[float] = None,
		inner: bool = False) -> pandas.DataFrame:
	"""
		Filters out timepoints that do not satisfy the detection criteria.
	Parameters
	----------
	left: pandas.Series
	right: pandas.Series
	dlimit: float
		Removes the timepoint if both points do no exceed this value.
	flimit: float
		If given, all points exceeding this value are treated as "undetected".
	inner: bool; default False
		If `true`, both points must exceed the detection limit for a given timepoint to be considered valid.

	Returns
	-------

	"""
	df = pandas.concat([left, right], axis = 1)
	if flimit is not None:
		# Mask the values so they are considered below the detection limit.
		# Assign a value of -1 so that they will be excluded wven if the detection limit is 0.
		left = left.mask(lambda s: s > flimit, -1)
		right = right.mask(lambda s: s > flimit, -1)
	df.columns = ['left', 'right']  # Prevent an IndexError when `left` and `right` refer to the same series.
	if inner:
		at_least_one_detected = (left > dlimit) & (right > dlimit)
	else:
		at_least_one_detected = (left > dlimit) | (right > dlimit)

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

		position_index_min = df.index.get_loc(position_index_min_value)
		position_index_max = df.index.get_loc(position_index_max_value)
		# Since we want to include the last index, increment position_index_max by one.
		position_index_max += 1
	result = df.iloc[position_index_min:position_index_max]
	return result


get_detected_points = get_valid_points


def format_linkage_matrix(Z, total_members: Optional[int]) -> pandas.DataFrame:
	linkage_dataframe = pandas.DataFrame(Z, columns = ["left", "right", "distance", "observations"])

	# linkage_dataframe.index = pandas.Index([i + len(squaremap.index) for i in linkage_dataframe.index], name = "clusterLabel")
	linkage_dataframe['left'] = linkage_dataframe['left'].astype(int)
	linkage_dataframe['right'] = linkage_dataframe['right'].astype(int)
	linkage_dataframe['observations'] = linkage_dataframe['observations'].astype(int)
	if total_members:
		linkage_dataframe['clusterId'] = list(total_members + i for i in range(len(linkage_dataframe)))
	return linkage_dataframe


def calculate_luminance(color: str) -> float:
	# 0.299 * color.R + 0.587 * color.G + 0.114 * color.B
	red = int(color[1:3], 16)
	green = int(color[3:5], 16)
	blue = int(color[5:], 16)

	lum = (.299 * red) + (.587 * green) + (.114 * blue)
	return lum / 255


def format_inconsistency_matrix(R) -> pandas.DataFrame:
	inconsistency_table = pandas.DataFrame(R, columns = ['mean', 'std', 'observations', 'statistic'])
	inconsistency_table['observations'] = inconsistency_table['observations'].astype(int)
	return inconsistency_table


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


if __name__ == "__main__":
	left = pandas.Series([0.00,0.00,0.00,0.00,0.00,0.263,0.07,0.081,0.069,0.042])
	right = pandas.Series([0.00,0.00,0.170,0.55,0.947,1.00,1.00,1.00,1.00,1.00])

	print(get_valid_points(left, right, 0.03, inner = True))
