from typing import Dict, List, Optional, Tuple, Mapping, Union

import pandas


class Cluster:
	""" Holds the possible ancestry candidates as well as the confidance score for each."""
	def __init__(self, initial_background: pandas.Series, timepoints: pandas.DataFrame):
		self.initial_background_label: str = initial_background.name
		self.timepoints = timepoints.copy()
		self.confidence: Dict[str,List[Tuple[str,float]]] = dict()
		initial_genotype = ['genotype-0']

		self.nests: Dict[str, List[str]] = {initial_background.name: initial_genotype}

	def add_genotype_to_background(self, unnested_label: str, nested_label: str, priority: Union[int, float]) -> None:
		if unnested_label not in self.nests:
			self.nests[unnested_label] = list()

		self.nests[unnested_label].append(nested_label)

		if unnested_label not in self.confidence:
			self.confidence[unnested_label] = []
		self.confidence[unnested_label].append((nested_label, priority))


	def get(self, label: str) -> List[str]:
		return self.nests[label]

	def is_a_member(self, label: str) -> bool:
		# return label in list(itertools.chain.from_iterable(self.nests.values()))
		return label in self.nests.keys()

	def get_sum_of_backgrounds(self) -> pandas.Series:
		background_labels = [k for k in self.nests.keys() if self.is_a_background(k)]
		background_frequencies = self.timepoints.loc[background_labels]
		# backgrounds = [v.trajectory for k, v in self.nests.items() if is_a_background(v)]
		# backgrounds = pandas.DataFrame(backgrounds)
		# total = backgrounds.sum()
		total = background_frequencies.sum()
		return total

	def is_a_background(self, element: str) -> bool:
		background = self.get(element)
		return len(background) == 1 or (len(background) == 2 and 'genotype-0' in background)
	def get_highest_priority(self, label:str)->Optional[str]:
		candidates = self.confidence.get(label, [])
		if candidates:
			# Explicity tell the sorting method to use the priority score.
			# This will prevent the method from using the genotype name to sort the elements,
			# So ties should be broken by whichever candidate was added as a candidate first.
			candidate, score = max(candidates, key = lambda s: s[1])
			if score < 1:
				candidate = None
		else:
			# `candidates` was an empty sequence.
			candidate = None
		return candidate


	def as_ancestry_table(self) -> pandas.Series:
		table = list()
		for identity, background in self.nests.items():
			# parent = self.get_highest_priority(identity)
			# if parent is None:
			parent = self.get_highest_priority(identity)
			if parent == identity or parent is None:
				parent = 'genotype-0'
			row = {
				'Parent':   parent,
				'Identity': identity
			}
			table.append(row)
		table = pandas.DataFrame(table)[['Parent', 'Identity']]  # Reorder columns

		return table.set_index('Identity')['Parent']

	def as_dict(self) -> Mapping[str, str]:
		return self.as_ancestry_table().to_dict()
