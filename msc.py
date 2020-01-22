from pathlib import Path
import random
import string
from tqdm import tqdm
if __name__ == "__main__":
	elements = list()

	alphabet = string.ascii_letters
	bar = tqdm(total = 10_000_000_000)
	for _ in range(0, 10_000_000_000):
		bar.update(1)
		l = random.choice(alphabet)
		r = random.choice(alphabet)
		elements.append((l, r))