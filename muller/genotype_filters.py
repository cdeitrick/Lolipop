import pandas

pandas.set_option('display.width', 300)
try:
	from muller import get_genotypes
	from muller.import_table import import_trajectory_table
except ModuleNotFoundError:
	# noinspection PyUnresolvedReferences
	import get_genotypes
	# noinspection PyUnresolvedReferences
	from import_table import import_trajectory_table
from pathlib import Path
from typing import Tuple


# noinspection PyTypeChecker
def get_invalid_genotype(genotypes: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float) -> pandas.DataFrame:
	fuzzy_fixed_cutoff = 0.5
	backgrounds = genotypes[genotypes.max(axis = 1) > fuzzy_fixed_cutoff]

	fuzzy_backgrounds = backgrounds.sum(axis = 1) > (1 + detection_cutoff)
	backgrounds = backgrounds[fuzzy_backgrounds]
	fixed_timepoints: pandas.DataFrame = backgrounds[backgrounds > fixed_cutoff].dropna(how = 'all').transpose()
	fixed_timepoints = fixed_timepoints.dropna(how = 'all').transpose()
	fixed_timepoints = [s.first_valid_index() for _, s in fixed_timepoints.iterrows()]
	not_backgrounds = genotypes[~genotypes.index.isin(backgrounds.index)]

	for background_id, background in backgrounds.iterrows():
		fuzzy_fixed_timepoints = background[background > fuzzy_fixed_cutoff]
		fixed_point: int = fuzzy_fixed_timepoints.first_valid_index()

		for genotype_label, genotype in not_backgrounds.iterrows():
			if genotype_label in backgrounds.index: continue

			detected_points = genotype[genotype > detection_cutoff]

			if not detected_points.empty:
				first_detected = detected_points.first_valid_index()
				last_detected = detected_points.last_valid_index()
			else:
				first_detected = 0
				last_detected = 0

			#print(first_detected, fixed_point, last_detected)
			if first_detected < fixed_point < last_detected:
				first_fixed_points: pandas.Series = genotype[fixed_timepoints]

				zero_points = first_fixed_points[first_fixed_points <= detection_cutoff]

				# The genotype was seen both before and after the background fixed.
				# Check if it was undetected when a genotype fixed.
				if zero_points.empty:
					return genotype_label


DF = pandas.DataFrame


def workflow(trajectories_filename: Path, detection_cutoff: float, fixed_cutoff: float = None) -> Tuple[DF, DF]:
	if fixed_cutoff is None: fixed_cutoff = 1 - detection_cutoff

	goptions = get_genotypes.GenotypeOptions.from_breakpoints(0.03)
	trajectory_table, _ = import_trajectory_table(trajectories_filename)
	genotype_table = get_genotypes.workflow(trajectories_filename,
		options = goptions)
	# return trajectory_table, genotype_table
	members = genotype_table.pop('members')

	for index in range(20):
		current_invalid_genotype = get_invalid_genotype(genotype_table, detection_cutoff, fixed_cutoff)
		if current_invalid_genotype is None:
			break
		else:
			invalid_members = members.loc[current_invalid_genotype].split('|')
			trajectory_table = trajectory_table[~trajectory_table.index.isin(invalid_members)]
			genotype_table = get_genotypes.workflow(trajectory_table, options = goptions)
			members = genotype_table.pop('members')
	else:
		print("Could not filter the genotypes after 20 iterations.")
	genotype_table['members'] = members
	return trajectory_table, genotype_table


if __name__ == "__main__":
	from pathlib import Path

	trajectories = Path("/home/cld100/Documents/github/muller_diagrams/Data files/B1_muller_try1.xlsx")
	workflow(trajectories, 0.03, 0.97)
