from typing import List, Optional


class PopulationGenotypes:
	def __init__(self, genotypes: List[List[str]] = None):
		if genotypes is None:
			genotypes = []
		self.genotypes = dict()
		for genotype in genotypes:
			self.create_genotype(genotype)

	def __len__(self) -> int:
		largest_index = max(self.genotypes.keys(), key = lambda s: int(s.split('-')[-1]))
		return int(largest_index.split('-')[-1])

	def to_list(self) -> List[List[str]]:
		return list(v for k, v in sorted(self.genotypes.items()))

	def merge_trajectories(self, left: str, right: str):
		genotype_left = self.get_genotype_from_trajectory(left)
		genotype_right = self.get_genotype_from_trajectory(right)

		if genotype_left and genotype_right:
			self.merge_genotypes(genotype_left, genotype_right)
		elif genotype_left:
			self.add_to_genotype(genotype_left, right)
		elif genotype_right:
			self.add_to_genotype(genotype_right, left)
		else:
			self.create_genotype([left, right])

	def add_to_genotype(self, genotype_label: str, new_member: str):
		self.genotypes[genotype_label].append(new_member)

	def create_genotype(self, members: List[str]):
		genotype_index = len(self.genotypes) + 1
		genotype_label = f"genotype-{genotype_index}"
		self.genotypes[genotype_label] = members

	def merge_genotypes(self, genotype_a: str, genotype_b: str):
		if genotype_a == genotype_b: return None
		genotype_a_members = self.genotypes[genotype_a]
		genotype_b_members = self.genotypes[genotype_b]

		new_genotype = genotype_a_members + genotype_b_members
		from pprint import pprint
		pprint(self.genotypes)
		self.create_genotype(new_genotype)

		self.genotypes.pop(genotype_a)
		self.genotypes.pop(genotype_b)

	def get_genotype_from_trajectory(self, trajectory_label: str) -> Optional[str]:
		for genotype_label, genotype_members in self.genotypes.items():
			if trajectory_label in genotype_members:
				return genotype_label
		return None
