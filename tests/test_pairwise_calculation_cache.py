import pytest
from muller_genotypes.metrics import PairwiseCalculation
import itertools

@pytest.fixture
def empty_cache():
	return PairwiseCalculation()

@pytest.fixture
def small_cache():
	elements = ['element1', 'element2', 'element3']
	combinations = itertools.combinations(elements, 2)
	pair_array = {k:'|'.join(k) for k in combinations}
	return PairwiseCalculation(pair_array)

def test_cache_empty(empty_cache, small_cache):
	assert bool(empty_cache) == False
	assert bool(small_cache) == True # For clarity.


