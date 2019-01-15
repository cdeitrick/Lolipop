import unittest
from inheritance.order import *

class TestOrderClusters(unittest.TestCase):
	def test_additie_check(self):
		left = pandas.Series([.3, .4, .1, .6])
		right = pandas.Series([.2, .4, .9, .4])

		result = check_additive_background(left, right, .05, .15)
		self.assertTrue(result)

if __name__ == "__main__":
	pass