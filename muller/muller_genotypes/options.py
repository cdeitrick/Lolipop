from typing import List, Optional

from dataclasses import dataclass


@dataclass
class GenotypeOptions:
	detection_breakpoint: float  # Minimum frequency to be considered detected.
	fixed_breakpoint: float  # Frequency at which a mutation is considered fixed.
	n_binom: Optional[int]
	similarity_breakpoint: float  # The cutoff indicating two trajectories are related.
	# The cutoff indicating two trajectories that were originally sorted into the same genotype are not
	# actually related.
	difference_breakpoint: float
	method: str

	@classmethod
	def from_matlab(cls) -> 'GenotypeOptions':
		return GenotypeOptions(
			detection_breakpoint = 0.03,
			fixed_breakpoint = 0.97,
			n_binom = 5,
			similarity_breakpoint = 0.05,
			difference_breakpoint = 0.10,
			method = 'matlab'
		)


@dataclass
class SortOptions:
	detection_breakpoint: float
	significant_breakpoint: float
	fixed_breakpoint: float
	frequency_breakpoints: List[float]  # Used when sorting non-fixed muller_genotypes.

	@classmethod
	def from_matlab(cls) -> 'SortOptions':
		return SortOptions(
			detection_breakpoint = 0.03,
			significant_breakpoint = 0.15,
			fixed_breakpoint = 0.85,
			frequency_breakpoints = [0.90, 0.75, 0.60, 0.45, 0.30, 0.15, 0.00]
		)
