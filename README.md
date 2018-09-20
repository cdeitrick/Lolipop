# A set of scripts to cluster mutations into genotypes and cluster genotypes by background.
![muller_plot](./example/B1_muller_try1.muller.png)

# Input Parameters

The script operates on a table listing all mutations and their corresponding frequencies at each timepoint or a table with each genotype and frequency at each timepoint.
The table must have a column named `Trajectory` or `Genotype` and integer columns for each timepoint. All other columns will be ignored. The frequencies can be represented as either a number between [0,1], a number between [0, 100] or as percentage.
The `Trajectory` and `Genotype` columns can contain any kind of label, but must be unique for each trajectory/genotype.

|            |            |            |          |       |          |   |       |       |       |       |       |       |
|------------|------------|------------|----------|-------|----------|---|-------|-------|-------|-------|-------|-------|
| Population | Trajectory | Chromosome | Position | Class | Mutation | 0 | 17    | 25    | 44    | 66    | 75    | 90    |
| B2         | 1          | 1          | 38102    | SNP   | C>T      | 0 | 0     | 26.1% | 100%  | 100%  | 100%  | 100%  |
| B2         | 2          | 1          | 62997    | SNP   | T>G      | 0 | 0     | 0     | 52.5% | 45.4% | 91.1% | 91%   |
| B2         | 3          | 1          | 78671    | SNP   | A>C      | 0 | 0     | 0     | 14.7% | 45%   | 92.4% | 88.7% |
| B2         | 4          | 1          | 96585    | SNP   | T>G      | 0 | 0     | 0     | 0     | 21.1% | 81.1% | 81.3% |
| B2         | 5          | 1          | 115010   | SNP   | G>T      | 0 | 0     | 0     | 40.3% | 48.9% | 5.7%  | 8%    |
| B2         | 6          | 1          | 156783   | SNP   | C>G      | 0 | 0     | 0     | 0     | 0     | 100%  | 100%  |
| B2         | 7          | 1          | 176231   | SNP   | T>A      | 0 | 0     | 0     | 27.3% | 78.1% | 100%  | 100%  |
| B2         | 8          | 1          | 205211   | SNP   | C>T      | 0 | 0     | 0     | 0     | 34.5% | 83.3% | 79.3% |
| B2         | 9          | 1          | 223199   | SNP   | C>G      | 0 | 0     | 0     | 0     | 0     | 26.9% | 34%   |
| B2         | 10         | 1          | 262747   | SNP   | T>C      | 0 | 0     | 11.7% | 0     | 0     | 0     | 10.3% |
| B2         | 11         | 1          | 264821   | SNP   | C>T      | 0 | 0     | 0     | 10.8% | 15.1% | 0     | 0     |
| B2         | 12         | 1          | 298548   | SNP   | G>A      | 0 | 12.5% | 0     | 15.3% | 18.1% | 17.5% | 19.1% |
| B2         | 13         | 1          | 299331   | SNP   | G>A      | 0 | 0     | 0     | 0     | 25.8% | 5.7%  | 7.5%  |
| B2         | 14         | 1          | 299332   | SNP   | C>T      | 0 | 38%   | 43.2% | 0     | 0     | 0     | 0     |
| B2         | 15         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 6.6%  | 10.4% | 6.2%  | 0     | 0     |
| B2         | 16         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 0     | 0     | 20.9% | 20.9% | 0     |
| B2         | 17         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 0     | 0     | 0     | 26.6% | 31.2% |
| B2         | 18         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 0     | 11.5% | 0     | 13.1% | 0     |
| B2         | 19         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 0     | 18.8% | 17.1% | 23.2% | 24.4% |
| B2         | 20         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 0     | 13.8% | 29.5% | 0     | 8.1%  |
| B2         | 21         | 1          | 299332   | SNP   | C>T      | 0 | 0     | 0     | 11.4% | 0     | 11%   | 12.3% |




# Sample Usage

```
python muller.py --input [input filename] --output [output folder]
```

The script generates the following files:
-  `input filename`.ggmuller_edges.tsv

	Used as the `edges` input to ggmuller.

- `input filename`.ggmuller_populations.tsv

	Used as the `population` input to ggmuller.

- `input filename`.genotypemembers.tsv

    A table indicating which trajectories compose each genotype.

- `input filename`.genotypes.original.tsv

    The initial genotypes generated before applying the genotype filters.

- `input filename`.genotypes.tsv

	A tab-delimited table with the mean frequency of each genotype at each timepoint. The mean is calculated from the trajectories that comprise each genotype.

- `input filename`.trajectories.original.tsv

    The unfiltered trajectories.

- `input filename`.trajectories.tsv

	A tab-delimited table of the population trajectories used in the analysis. Each trajectory represents the frequency of a single mutation at each timepoint.

- `input filename`.mermaid.md

	A script written in the [mermaid](https://mermaidjs.github.io) scripting language. Generates a diagram indicating the hierarchy of genotypes/backgrounds in the current population.

- `input filename`.mermaid.png

    If [mermaid.cli](https://github.com/mermaidjs/mermaid.cli) is installed, the mermaid script will automatically be used to generate a map of nested genotypes.

- `input filename`.yaml

	Contains additional information for each genotype, as well as the parameters used for the analysis.

- `input filename`.r

	A basic r script to import the `population` and `edges` tables and generate a muller plot.
	
- `input filename`.png

	The muller plot generated by the r script. Colors correspond to the same genotypes in the genotype plots.

- `input_filename`.genotypeplot.png

    Plots of all trajectories (colored by parent genotype) and genotypes.


# Script Options

	-h, --help                  Show this help message and exit
	-i, --input                 The table of trajectories to cluster. Must be an excel file or csv/tsv file.
	-o,  --output               The folder to save the files to.
	--genotypes                 Indicates that the input table contains genotypes rather
	                            than mutational trajectories.
	-u, --uncertainty           The uncertainty to apply when performing
	                            frequency-based calculations. For
	                            example, a frequency at a given timepoint
	                            is considered undetected if it falls
	                            below 0 + `uncertainty`.
	--fixed                     The minimum frequency at which to
	                            consider a mutation fixed. Defaults to 
	                            1 - `uncertainty`
	-s, --significant           [0.15] The frequency at which to consider a genotype 
	                            significantly greater than zero.
	--matlab                    Mimics the output of the original matlab script.
	-f, --frequencies           [0.10] The frequency cutoff to use when sorting genotypes. 
	                            May be a comma-separated string of frequencies, or a set inverval 
	                            to use when generating the frequency breakpoints. 
	                            For example, a value of 0.15 will use the frequencies 0,.15,.30,.45...
	-r --similarity-cutoff      [0.05] Maximum p-value difference to consider trajectories related. 
	                            Used when grouping trajectories into genotypes.
	-l, --difference-cutoff     [0.10] Used to unlink unrelated trajectories present in a genotype.

# Example Usage
```
python --input [filename] --output [folder] --matlab
```
Trajectories will be grouped into genotypes and each genotype will be nested using the same parameters the original matlab script used. The original script used different values in each script when performing frequency detection/significance. The improved version of the script (without the --matlab flag) harmonized these parameters so that similar calculations use the same parameter rather than re-defining the value to use.
```
python --input [filename] --frequencies 0.05 --detected 0.10
```
Groups genotypes in groups of 0.05 (i.e. [0.00, 0.05, 0.10, ... , 0.90, 0.95, 1.00]) based on each genotype's maximum frequency. Each genotype in each group is then sorted by the timepoint it was first detected (the first timepoint where the frequency was greater than 0.10). Output files are saved to the same folder as the input table.

# Genotype Plots
The `.genotypeplot.png` file gives an easy-to-visualize plot of all trajectory and genotype frequencies over time.
Trajectories are colored based on their parent genotype.
![genotypeplot](./example/example.genotypeplot.png)

# Mermaid Diagram

The `.mermaid` file can be used to generate a quick diagram showing the relation between all genotypes in the population.

![diagram](./example/sample_mermaid_diagram.png)


