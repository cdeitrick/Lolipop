import pandas
pandas.set_option('display.width', 300)
try:
	from muller import get_genotypes
except ModuleNotFoundError:
	import get_genotypes

def filter_genotypes(trajectories:pandas.DataFrame, fixed_cutoff:float, detection_cutoff:float)->pandas.DataFrame:

	genotypes = get_genotypes.workflow(trajectories, options = get_genotypes.GenotypeOptions.from_matlab())
	members = genotypes.pop('members')
	backgrounds = genotypes[genotypes.max(axis = 1) > fixed_cutoff]

	for genotype_label, genotype in genotypes.iterrows():
		for background_label, background in backgrounds.iterrows():

			difference = background - genotype

			genotype_above_background = difference[difference > detection_cutoff]






if __name__ == "__main__":
	from pathlib import Path
	trajectories = Path("/home/cld100/Documents/github/muller_diagrams/Data files/B1_muller_try1.xlsx")
	filter_genotypes(trajectories, 0.97, 0.03)
