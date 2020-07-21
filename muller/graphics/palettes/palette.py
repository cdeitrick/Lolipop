from pathlib import Path
from typing import Dict, List, Optional


class Palette:
	# Holds the palette used by genotypes/trajectories.
	def __init__(self, name:str, palette: Dict[str, str], members: Optional[Dict[str, List[str]]] = None):
		self.name = name
		self.palette_genotype = palette
		if members:
			self._palette_trajectory = self.generate_trajectory_palette(palette, members)
		else:
			self._palette_trajectory = None

	def __call__(self, item:str, default = None)->str:
		if item in self.palette_genotype:
			return self.palette_genotype[item]
		elif self._palette_trajectory is not None and item in self._palette_trajectory:
			return self._palette_trajectory[item]
		else:
			return default

	def __getitem__(self, item):
		return self(item)

	def generate_trajectory_palette(self, palette: Dict[str, str], members: Dict[str, List[str]]) -> Dict[str, str]:
		""" Maps the genotype palette to each of the member trajectories."""
		p = dict()
		for genotype_name, members in members.items():
			for member in members:
				p[member] = palette[genotype_name]
		self._palette_trajectory = p
		return p

	def get_genotype_palette(self) -> Dict[str, str]:
		return self.palette_genotype

	def get_trajectory_palette(self) -> Dict[str, str]:
		if self._palette_trajectory is None:
			message = f"The trajectory palette was not generated."
			raise ValueError(message)
		return self._palette_trajectory

	def get(self, item:str, default = None)->str:
		return self.palette_genotype.get(item, self._palette_trajectory.get(item, default))

	def save(self, filename:Path):
		data = {
			'name': self.name,
			'genotype': self.palette_genotype,
			'trajectory': self._palette_trajectory
		}
		import json
		filename.write_text(json.dumps(data, indent = 4, sort_keys = True))

if __name__ == "__main__":
	pass
