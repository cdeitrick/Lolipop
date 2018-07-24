
if __name__ == "__main__":
	left = [1,2,3]
	right = [4,5,6]

	center = [
		[1,5,3],
		[1,2,3],
		[1,2,3,4,5,6],
		[4,5,6]
	]

	aa = center[0]
	aa.append(897)
	center.remove(right)
	from pprint import pprint

	pprint(center)