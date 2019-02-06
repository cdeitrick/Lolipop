import pandas
from dataclasses import dataclass
from typing import List, Dict
import itertools
@dataclass
class Genotype:
	name: str
	background: List[str]


class Cluster:
	def __init__(self, initial_background: pandas.Series, timepoints:pandas.DataFrame):
		self.initial_background_label: str = initial_background.name
		self.timepoints = timepoints.copy()
		initial_genotype = Genotype(
			name = initial_background.name,
			background = [initial_background.name]
		)
		initial_genotype = ['genotype-0']

		self.nests:Dict[str,List[str]] = {initial_background.name: initial_genotype}

	def add_genotype_to_background(self, unnested_label:str, nested_label: str)->None:
		if unnested_label in self.nests:
			self.nests[unnested_label].append(nested_label)
		else:
			self.nests[unnested_label] = [nested_label]
	def get(self, label: str) -> Genotype:
		return self.nests[label]

	def is_a_member(self, label: str) -> bool:
		#return label in list(itertools.chain.from_iterable(self.nests.values()))
		return label in self.nests.keys()

	def get_sum_of_backgrounds(self) -> pandas.Series:
		background_labels = [k for k in self.nests.keys() if self.is_a_background(k)]
		background_frequencies = self.timepoints.loc[background_labels]
		#backgrounds = [v.trajectory for k, v in self.nests.items() if is_a_background(v)]
		#backgrounds = pandas.DataFrame(backgrounds)
		#total = backgrounds.sum()
		total = background_frequencies.sum()
		return total

	def is_a_background(self, element: str) -> bool:
		background = self.get(element)
		return len(background) == 1 or (len(background) == 2 and 'genotype-0' in background)
