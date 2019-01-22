import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas

PairwiseArrayType = Dict[Tuple[str, str], float]


class PairwiseCalculation:
	"""
		Calculates and holds the calculations for all pairwise elements.

		Usage
		-----
		calculation = PairwiseCalculation()
		# Update values
		pair_array = calculation.update_values(timeseries, detection_cutoff, fixed_cutoff)
	"""

	def __init__(self, pairwise_array: PairwiseArrayType = None):
		if pairwise_array is None:
			self.pairwise_values = dict()
		else:
			self.pairwise_values: PairwiseArrayType = pairwise_array

	def __bool__(self) -> bool:
		return bool(self.pairwise_values)

	def asdict(self) -> PairwiseArrayType:
		return self.pairwise_values

	def squareform(self):
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

	def get(self, left, right, default = None) -> float:
		result = self.pairwise_values.get((left, right), default)
		return result

	def reduce(self, labels: Iterable[str]) -> 'PairwiseCalculation':
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

	def update(self, pair_array: PairwiseArrayType) -> 'PairwiseCalculation':
		self.pairwise_values.update(pair_array)
		return self

	def unique(self) -> List[Tuple[str, str]]:
		seen = set()
		for key in self.pairwise_values.keys():
			if key in seen or key[::-1] in seen:
				continue
			else:
				yield key

	def save(self, filename: Optional[Path]):
		with filename.open('w') as output:
			for key in self.unique():
				value = self.get(*key)
				line = f"{key[0]}\t{key[1]}\t{value}\n"
				output.write(line)

	@classmethod
	def read(cls, filename: Path) -> 'PairwiseCalculation':
		contents = filename.read_text().split('\n')
		contents = [i.split('\t') for i in contents]
		data = dict()
		for line in contents:
			left, right, value = line
			value = float(value)
			data[left, right] = value
			data[right, left] = value
		return PairwiseCalculation(data)
