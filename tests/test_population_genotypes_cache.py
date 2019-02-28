import pytest

from clustering.methods.population_genotypes_cache import PopulationGenotypes


@pytest.fixture
def empty_genotypes() -> PopulationGenotypes:
	return PopulationGenotypes()


@pytest.fixture
def one_genotype() -> PopulationGenotypes:
	return PopulationGenotypes([["2", "3"]])


@pytest.fixture
def many_genotypes() -> PopulationGenotypes:
	g = [
		["1", "4", "5"],
		["2"],
		["7", "9"],
		["10"],
		["12"]
	]
	return PopulationGenotypes(g)


def test_create_genotypes(empty_genotypes):
	members = ["1", "4"]

	empty_genotypes.create_genotype(members)
	expected = {
		'genotype-1': members
	}
	assert expected == empty_genotypes.genotypes

	empty_genotypes.create_genotype(["5"])
	expected = {
		'genotype-1': members,
		'genotype-2': ["5"]
	}
	assert expected == empty_genotypes.genotypes


def test_get_genotype_from_trajectory(many_genotypes):
	assert 'genotype-3' == many_genotypes.get_genotype_from_trajectory("9")
	assert 'genotype-5' == many_genotypes.get_genotype_from_trajectory("12")
	assert None == many_genotypes.get_genotype_from_trajectory(12)


def test_add_to_genotype(many_genotypes):
	many_genotypes.add_to_genotype("genotype-1", "15")

	assert many_genotypes.genotypes['genotype-1'] == ["1", "4", "5", "15"]


def test_length(many_genotypes):
	many_genotypes.genotypes['genotype-15'] = ["asd"]
	assert len(many_genotypes) == 15


def test_merge_genotypes(many_genotypes):
	expected = {
		'genotype-1': ["1", "4", "5"],
		'genotype-3': ["7", "9"],
		'genotype-4': ["10"],
		'genotype-6': ["2", "12"]
	}
	many_genotypes.merge_genotypes("genotype-2", "genotype-5")
	assert expected == many_genotypes.genotypes


def test_merge_trajectories(one_genotype):
	one_genotype.merge_trajectories("2", "1")
	expected = {'genotype-1': ["2", "3", "1"]}
	assert expected == one_genotype.genotypes

	one_genotype.merge_trajectories("17", "18")
	expected = {'genotype-1': ["2", "3", "1"], "genotype-2": ["17", "18"]}
	assert expected == one_genotype.genotypes

	one_genotype.merge_trajectories("17", "1")
	expected = {'genotype-3': ['17', '18', '2', '3', '1']}
	assert expected == one_genotype.genotypes

def test_to_list(many_genotypes):
	expected = [
		["1", "4", "5"],
		["2"],
		["7", "9"],
		["10"],
		["12"]
	]
	assert expected == many_genotypes.to_list()