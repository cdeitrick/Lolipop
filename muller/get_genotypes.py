
import pandas

import itertools
from pprint import pprint
import math
from dataclasses import dataclass
from typing import List, Optional
from pathlib import Path
try:
	from muller.time_series_import import import_timeseries
except ModuleNotFoundError:
	from time_series_import import import_timeseries

@dataclass
class PairArrayValue:
	"""
		Holds the values assigned to parray in the original script.
	"""
	population: int
	left: int
	right: int
	p_value: float
	sigma_value: float
	difbar: float


def calculate_pairwise_similarity(population_id: int, timeseries: pandas.DataFrame) -> List[PairArrayValue]:
	trajectories = timeseries[[i for i in timeseries.columns if i not in ['Population', 'Position', 'Trajectory']]]
	combos = sorted(itertools.combinations(timeseries['Trajectory'].values, 2))
	pair_array = list()
	for pair in combos:
		left, right = pair
		# trajectory indicies start at 1, need to convert to 0-indexed.
		left_trajectories = trajectories.iloc[left - 1]
		right_trajectories = trajectories.iloc[right - 1]

		df = pandas.concat([left_trajectories, right_trajectories], axis = 1)

		filtered_df = df[(df < 0.97).any(axis = 1)]
		filtered_df = filtered_df[(filtered_df > 0.03).any(axis = 1)]
		if filtered_df.empty:
			pair_array_value = PairArrayValue(
				population_id, pair[0], pair[1], 0.0, math.nan, math.nan
			)
		else:
			ps = filtered_df.mean(axis = 1)
			sigma2 = (ps * (1 - ps)) / 5
			difference = filtered_df.iloc[:, 0] - filtered_df.iloc[:, 1]
			sigmapair = math.sqrt(sigma2.sum()) / len(difference)
			difbar = (abs(difference) / len(difference)).sum()
			pval = 1 - math.erf(difbar / (math.sqrt(2) * sigmapair))

			pair_array_value = PairArrayValue(
				population_id, pair[0], pair[1], pval, sigmapair, difbar
			)
		pair_array.append(pair_array_value)
	return pair_array


def find_genotype(element: int, genotypes: List[List[int]]) -> List[int]:
	candidates = [i for i in genotypes if element in i]
	try:
		value = candidates[0]
	except IndexError:
		value = None

	return value


def find_all_genotypes(pairs: List[PairArrayValue], relative_cutoff: float = 0.05) -> List[List[int]]:
	genotypes = [[1]]
	for pair in pairs:
		if pair.p_value > relative_cutoff:  # are the genotypes related?
			# Check if any of the trajectories are already listed in genotypes.
			genotype_left = find_genotype(pair.left, genotypes)
			genotype_right = find_genotype(pair.right, genotypes)

			if genotype_left and genotype_right:
				if genotype_left != genotype_right:  # they are listed under two different genotypes. Combine them.
					genotype_left += genotype_right
					# Remove the redundant genotypes.
					genotypes.remove(genotype_right)
			elif genotype_left:
				genotype_left.append(pair.right)
			elif genotype_right:
				genotype_right.append(pair.left)
			else:  # Neither element is listed. Create a new genotype
				genotypes.append([pair.left, pair.right])
	return genotypes


def get_p_value(left: int, right: int, pairwise_array: List[PairArrayValue]) -> Optional[PairArrayValue]:
	for candidate in pairwise_array:
		l = candidate.left == left and candidate.right == right
		r = candidate.left == right and candidate.right == left
		if l or r:
			return candidate
	else:
		return None


def split_genotype_in_two(genotype: List[int], unlinked_trajectories: pandas.DataFrame,
		pair_array: List[PairArrayValue], link_cut: float):
	minimum_pvalue_index = (unlinked_trajectories['pvalue'] - link_cut).abs().idxmin()
	minimum_pvalue_genotype = unlinked_trajectories.loc[minimum_pvalue_index]
	#print(minimum_pvalue_genotype)
	new_genotype_1_base = int(minimum_pvalue_genotype['left'])
	new_genotype_2_base = int(minimum_pvalue_genotype['right'])
	new_genotype_1_base, new_genotype_2_base = sorted([new_genotype_1_base, new_genotype_2_base])
	new_genotype_1 = [new_genotype_1_base]
	new_genotype_2 = [new_genotype_2_base]

	for genotype_member in genotype:
		genotype_member = int(genotype_member)
		# if genotype_member in new_genotype_1 or genotype_member in new_genotype_2:
		if genotype_member in new_genotype_1 or genotype_member in new_genotype_2:
			pass
		else:
			p_value_1 = get_p_value(new_genotype_1_base, genotype_member, pair_array)
			p_value_2 = get_p_value(new_genotype_2_base, genotype_member, pair_array)

			if p_value_1.p_value >= p_value_2.p_value:
				new_genotype_1.append(genotype_member)
			else:
				new_genotype_2.append(genotype_member)
	return new_genotype_1, new_genotype_2


def split_unlinked_genotypes(genotypes: List[List[int]], pair_array: List[PairArrayValue], link_cutoff: float) -> List[
	List[int]]:
	for genotype in genotypes:  # Operate on a copy of genotypes
		if len(genotype) > 1:
			combination_pairs = list()
			for combination_pair in itertools.combinations(genotype, 2):
				left, right = combination_pair
				p_value = get_p_value(left, right, pair_array)
				value = [left, right, p_value.p_value]
				combination_pairs.append(value)
			genotype_combinations = pandas.DataFrame(combination_pairs, columns = ['left', 'right', 'pvalue'])
			#print(genotype_combinations)
			unlinked_trajectories = genotype_combinations[genotype_combinations['pvalue'] <= link_cutoff]

			if len(unlinked_trajectories) != 0:
				new_genotype_1, new_genotype_2 = split_genotype_in_two(genotype[:], genotype_combinations, pair_array,
					link_cutoff)

				genotypes.append(new_genotype_1)
				genotypes.append(new_genotype_2)
				genotypes.remove(genotype)
	return sorted(genotypes, key = lambda s: len(s))


def get_genotypes(timeseries: pandas.DataFrame):
	npops = 1
	plotting = 0
	relative_cutoff = 0.10
	link_cutoff = 0.05
	totals = 0

	#colors = varycolor(len(timeseries))
	populations = timeseries.groupby(by = 'Population')
	all_genotypes = dict()
	for population_id, population_data in populations:
		# Trajectories represent the population frequencies at each timepoint
		# Each row represents a single timepoint, each column represents a mutation.

		pair_array = calculate_pairwise_similarity(population_id, timeseries)

		genotypes = find_all_genotypes(pair_array, relative_cutoff)

		# at the end, look at all trajectories that are not listed and
		# append them as their own category.
		flattened_genotypes = list(itertools.chain.from_iterable(genotypes))
		other_trajectories = [i for i in timeseries['Trajectory'].values if i not in flattened_genotypes]
		genotypes.append(other_trajectories)

		# finally, for each genotype, make sure each trajectory pair has some
		# non-trivial linkage (say, >0.0005). this avoids falsely linking together
		# trajectories. if not, divide the offending trajectories into two
		# camps, and sort the rest according to which one they are more
		# closely linked to. repeat until everything is linked.
		# pprint(pair_array)

		while True:
			starting_size_of_the_genotype_array = len(genotypes)
			genotypes = split_unlinked_genotypes(genotypes[:], pair_array, link_cutoff)
			if len(genotypes) == starting_size_of_the_genotype_array:
				break

		all_genotypes[population_id] = genotypes
	# TODO Fix so that it returns a genotype for each population
	return genotypes
def get_mean_genotypes(genotypes, timeseries):
	mean_genotypes = list()
	for genotype in genotypes:
		genotype_timeseries = timeseries[timeseries['Trajectory'].isin(genotype)]
		mean_genotype_timeseries = genotype_timeseries.mean()
		mean_genotype_timeseries.name = "|".join(map(str,genotype))
		mean_genotypes.append(mean_genotype_timeseries)

	mean_genotypes = pandas.DataFrame(mean_genotypes)
	if 'Trajectory' in mean_genotypes:
		mean_genotypes.pop('Trajectory')
	if 'Position' in mean_genotypes:
		mean_genotypes.pop('Position')
	return mean_genotypes

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
	'-i', '--input',
	help = "The table of trajectories to cluster.",
	action = 'store',
	dest = 'filename'
)
@dataclass
class Arguments:
	filename: str
	@classmethod
	def from_parser(cls, parser):
		return Arguments(parser.filename)
if __name__ == "__main__":
	from muller.variables import filename
	#parser = parser.parse_args()
	parser = Arguments(filename)
	filename = Path(parser.filename)
	timepoints, info = import_timeseries(filename)

	genotypes = get_genotypes(timepoints)
	mean_genotypes = get_mean_genotypes(genotypes, timepoints)

	mean_genotypes.to_csv(str(filename.with_suffix('.mean.tsv')), sep = '\t')


