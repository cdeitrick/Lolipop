from typing import Dict, List, Optional
import logging
logger = logging.getLogger(__file__)
import pandas



class Cluster:
	def __init__(self, initial_background: pandas.Series, timepoints: pandas.DataFrame):
		self.initial_background_label: str = initial_background.name
		self.timepoints = timepoints.copy()
		self.confidence = dict()
		initial_genotype = ['genotype-0']

		self.nests: Dict[str, List[str]] = {initial_background.name: initial_genotype}

	def add_genotype_to_background(self, unnested_label: str, nested_label: str, priority:int = None) -> None:
		logger.info(f"adding {nested_label} as a potential background for {unnested_label} with priority {priority}")
		if unnested_label not in self.nests:
			self.nests[unnested_label] = list()

		if priority and False:
			self.nests[unnested_label].insert(0, nested_label)
		else:
			self.nests[unnested_label].append(nested_label)
		self.nests[unnested_label].append(nested_label)

		if priority:
			self.confidence[nested_label, unnested_label] = priority
			# To make lookup easier
			self.confidence[unnested_label, nested_label] = priority

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
		candidates = dict()
		for (left, right), value in self.confidence.items():
			if left == label:
				candidates[right] = value
			elif right == label:
				candidates[left] = value
		try:
			maximum = max(candidates.values())
		except ValueError:
			maximum = None
		try:
			candidate_backgrounds = [i for i in self.get(label) if candidates.get(i) == maximum]
			return candidate_backgrounds[0]
		except IndexError:
			return None

	def as_ancestry_table(self) -> pandas.Series:
		table = list()
		for identity, background in self.nests.items():
			#parent = self.get_highest_priority(identity)
			#if parent is None:
			parent = background[0]
			if parent == identity:
				parent = 'genotype-0'

			row = {
				'Parent':   parent,
				'Identity': identity
			}
			table.append(row)
		table = pandas.DataFrame(table)[['Parent', 'Identity']]  # Reorder columns

		return table.set_index('Identity')['Parent']
	def as_dict(self)->Dict[str,str]:
		return self.as_ancestry_table().to_dict()