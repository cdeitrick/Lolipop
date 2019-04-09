import pytest
import treetools


@pytest.mark.parametrize(
	"left,right,expected",
	[
		("abcdef","defjkl", 3),
		("abcdef", "zyx", 0),
		("jklmnop", "mno", 3)
	]
)
def test_common_values(left, right, expected):
	assert treetools.common_elements(left, right) == expected

def test_tokenize():
	expected = ["abc", "def", "jkl", "mno", "pqr"]

	result = treetools.tokenize(["abc def", "jkl", "mno pqr"])

	assert result == expected
