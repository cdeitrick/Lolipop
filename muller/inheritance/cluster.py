import pandas
from dataclasses import dataclass
from typing import List
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

		self.nests = {initial_background.name: initial_genotype}

	def add_genotype_to_background(self, genotype_label: str, type_genotype: Genotype) -> None:
		if genotype_label in self.nests:
			self.nests[genotype_label].background += type_genotype.background

			if len(self.nests[genotype_label].background) > 2 and self.initial_background_label in self.nests[genotype_label].background:
				self.nests[genotype_label].background.remove(self.initial_background_label)

		else:
			self.nests[genotype_label] = type_genotype

	def get(self, label: str) -> Genotype:
		return self.nests[label]

	def is_a_member(self, label: str) -> bool:
		return label in itertools.chain.from_iterable(i.background for i in self.nests.values())

	def get_sum_of_backgrounds(self) -> pandas.Series:
		background_labels = [v.name for v in self.nests.values() if is_a_background(v)]
		background_frequencies = self.timepoints.loc[background_labels]
		#backgrounds = [v.trajectory for k, v in self.nests.items() if is_a_background(v)]
		#backgrounds = pandas.DataFrame(backgrounds)
		#total = backgrounds.sum()
		total = background_frequencies.sum()
		return total

def is_a_background(element: Genotype) -> bool:
	return len(element.background) == 1
