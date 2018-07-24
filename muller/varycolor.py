
import random
from typing import List
def generate_random_color()->str:
	r = random.randint(100,255)
	g = random.randint(100,255)
	b = random.randint(100,255)

	color = "#{:>02X}{:>02X}{:>02X}".format(r, g, b)
	return color

def varycolor(number_of_plots:int)->List[str]:
	return [generate_random_color() for _ in range(number_of_plots)]