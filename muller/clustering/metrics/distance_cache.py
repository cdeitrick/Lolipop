import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas
from scipy.spatial import distance

PairwiseArrayType = Dict[Tuple[str, str], float]

# TODO: Refactor using UserDict
class DistanceCache:
	"""
		Calculates and holds the distances for all pairwise elements.

		Usage
		-----
		calculation = PairwiseCalculation()
		# Update values
		pair_array = calculation.update_values(timeseries, detection_cutoff, fixed_cutoff)
	"""

	def __init__(self, pairwise_array: PairwiseArrayType = None):
		if pairwise_array is None: pairwise_array = dict()
		self.pairwise_values = dict()

		for (left, right), value in pairwise_array.items():
			self.pairwise_values[left, right] = value
			self.pairwise_values[right, left] = value

	def __bool__(self) -> bool:
		return bool(self.pairwise_values)

	def __len__(self) -> int:
		return len(self.pairwise_values)
	def __getitem__(self, item):
		return self.pairwise_values[item]
	def asdict(self) -> PairwiseArrayType:
		return self.pairwise_values

	def squareform(self)->pandas.DataFrame:
		""" Converts a dictionary with all pairwise values for a set of points into a square matrix representation.
		"""
		keys = sorted(set(itertools.chain.from_iterable(self.pairwise_values.keys())))
		_square_map = dict()
		for left in keys:
			series = dict()
			for right in keys:
				value = self.get(left, right, 0)
				series[right] = value
			_square_map[left] = series
		return pandas.DataFrame(_square_map)

	def triangle(self):
		""" Returns the condensed squareform of the pair array."""
		# noinspection PyTypeChecker
		return distance.squareform(self.squareform().values)

	def get(self, left, right, default = None) -> float:
		try:
			result = self.pairwise_values[left, right]
		except KeyError:
			result = self.pairwise_values.get((right, left), default)
		return result

	def reduce(self, labels: Iterable[str]) -> 'DistanceCache':
		"""
			Removes all labels from `self.pairwise_values` that are not present in `labels`
		"""
		labels = set(labels)
		pair_array = {
			k: v for k, v in self.pairwise_values.items()
			if (k[0] in labels and k[1] in labels)
		}
		self.pairwise_values = pair_array
		return self

	def update(self, pair_array: PairwiseArrayType) -> 'DistanceCache':
		for (left, right), value in pair_array.items():
			self.pairwise_values[left, right] = value
			self.pairwise_values[right, left] = value
		return self

	def unique(self) -> List[Tuple[str, str]]:
		seen = set()
		for key in sorted(self.pairwise_values.keys()):
			if key in seen or key[::-1] in seen:
				continue
			else:
				seen.add(key)
				yield key

	def save(self, filename: Path):
		with filename.open('w') as output:
			for key in self.unique():
				value = self.get(*key)
				line = f"{key[0]}\t{key[1]}\t{value}\n"
				output.write(line)

	@property
	def values(self)->List[float]:
		return list(self.pairwise_values.values())

	@classmethod
	def read(cls, filename: Path) -> 'DistanceCache':
		contents = filename.read_text().split('\n')
		contents = [i.split('\t') for i in contents]
		data = dict()
		for line in contents:
			left, right, value = line
			value = float(value)
			data[left, right] = value
			data[right, left] = value
		return DistanceCache(data)

	@classmethod
	def from_squareform(cls, square:pandas.DataFrame)->'DistanceCache':
		data = dict()

		for left, row in square.iterrows():
			for right, value in row.items():
				data[left,right] = value
				data[right,left] = value
		return DistanceCache(data)
