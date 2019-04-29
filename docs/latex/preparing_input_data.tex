 # Preparing the Data
[Breseq](http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing) is a variant caller used to analyze samples genomes during microbial evolution experiments. An evolution experiment which sequences a population at selected timepoints will end up with a number of output tables reporting detected mutations and frequency of each mutation at each timepoint. Each table will look similar to this:


| evidence | position | mutation | freq  | annotation         | gene          | description                      |
|----------|----------|----------|-------|--------------------|---------------|----------------------------------|
| RA       | 14,350   | (A)9>8   | 1.30% | intergenic(48/298) | pepX_1</>dnaE | DNA polymerase III subunit alpha |
| RA       | 14,350   | (A)9>10  | 0.90% | intergenic(48/298) | pepX_1</>dnaE | DNA polymerase III subunit alpha |
| RA       | 19,123   | G>A      | 1.10% | E76K(GAA>AAA)      | pyk           | Pyruvate kinase                  |
| RA       | 23,912   | (A)8>9   | 1.00% | coding(243/312nt)  | PROKKA_00020  | hypothetical protein             |


To prepare the data for use by the muller scripts, add a new column to each table indicating the timepoint it represents, then combine all of the tables into a single table. The result should look something like this:

| timepoint | evidence | position | mutation | freq  | annotation         | gene          | description                      |
|-----------|----------|----------|----------|-------|--------------------|---------------|----------------------------------|
| 10        | RA       | 14,350   | (A)9>8   | 1.30% | intergenic(48/298) | pepX_1</>dnaE | DNA polymerase III subunit alpha |
| 10        | RA       | 14,350   | (A)9>10  | 0.90% | intergenic(48/298) | pepX_1</>dnaE | DNA polymerase III subunit alpha |
| 10        | RA       | 19,123   | G>A      | 1.10% | E76K(GAA>AAA)      | pyk           | Pyruvate kinase                  |
| 3         | RA       | 22,606   | C>A      | 1.50% | Q148K(CAG>AAG)     | PROKKA_00018  | hypothetical protein             |
| 9         | RA       | 23,912   | (A)8>9   | 2.50% | coding(243/312nt)  | PROKKA_00020  | hypothetical protein             |
| 10        | RA       | 23,912   | (A)8>9   | 1.00% | coding(243/312nt)  | PROKKA_00020  | hypothetical protein             |
| 6         | RA       | 26,051   | C>A      | 1.80% | S166I(AGT>ATT)     | PROKKA_00023< | putative permease                |
| 2         | RA       | 26,291   | (T)7>6   | 1.60% | coding(257/906nt)  | PROKKA_00023< | putative permease                |
| 2         | RA       | 35,252   | T>A      | 1.00% | E26V(GAA>GTA)      | PROKKA_00033< | Integrase core domain protein    |


Then, the data should be converted into a table such that each row represents a single mutation, and each timepoint occupies a unique column. An additional column, `Trajectory`, should then be added with a unique identifier for each mutation. The specific format of this identifier is arbitrary, and it may be convienient to simply number the mutations. It does not matter whether additional columns are included in the table (although they may be used to add for annotations later), as long as there is a `Trajectory` column with the unique identifiers for each mutation along with numeric columns for each sampled timepoint.

| Trajectory | position | mutation | gene                       | annotation          | description                                                        | 0    | 1    | 2    | 3    | 4    | 5    | 6    | 7    | 8    | 9    | 10   |
|------------|----------|----------|----------------------------|---------------------|--------------------------------------------------------------------|------|------|------|------|------|------|------|------|------|------|------|
| 1          | 36,414   | G>A      | speA                       | C127Y(TGT>TAT)      | Arginine decarboxylase                                             | 0.00 | 0.06 | 0.10 | 0.17 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 2          | 138,043  | C>T      | rlmCD_1                    | H441H(CAC>CAT)      | 23S rRNA (uracilC(5))methyltransferase RlmCD                       | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.02 | 0.05 | 0.02 |
| 3          | 165,470  | C>A      | PROKKA_00173/>PROKKA_00174 | intergenic(+174/91) | Relaxase/Mobilisation nuclease domain protein/hypothetical protein | 0.09 | 0.11 | 0.16 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 4          | 234,888  | C>A      | gapN                       | A95E(GCA>GAA)       |Some definitions NADPdependent glyceraldehyde3phosphate dehydrogenase               | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.25 | 0.39 | 0.22 | 0.40 | 0.14 | 0.03 |
| 5          | 264,552  | T>C      | rbgA                       | D241D(GAT>GAC)      | Ribosome biogenesis GTPase A                                       | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.03 | 0.09 | 0.04 | 0.02 | 0.00 |
| 6          | 407,633  | G>T      | gdhA                       | L206F(TTG>TTT)      | NADPspecific glutamate dehydrogenase                               | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.08 | 0.02 | 0.24 | 0.31 |
| 7          | 458,680  | A>C      | rlmI<                      | L58V(TTG>GTG)       | Ribosomal RNA large subunit methyltransferase I                    | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.08 | 0.09 | 0.07 |
| 8          | 636,386  | G>A      | lpd_2                      | R34H(CGT>CAT)       | Dihydrolipoyl dehydrogenase                                        | 0.00 | 0.00 | 0.00 | 0.00 | 0.15 | 0.38 | 0.34 | 0.29 | 0.33 | 0.15 | 0.03 |
| 9          | 693,913  | C>A      | PROKKA_00705<              | A116S(GCC>TCC)      | putative ABC transporter ATPbinding protein/MT1014                 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.06 | 0.06 | 0.11 |
| 10         | 701,443  | T>G      | ileS<                      | K107T(AAG>ACG)      | IsoleucinetRNA ligase                                              | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.06 | 0.07 | 0.04 |
| 11         | 764,481  | G>A      | bglF_1<                    | I192I(ATC>ATT)      | PTS system betaglucosidespecific EIIBCA component                  | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.02 | 0.02 | 0.04 | 0.05 | 0.03 |
| 12         | 860,048  | C>A      | arnB<                      | R288L(CGC>CTC)      | UDP4amino4deoxyLarabinoseoxoglutarate aminotransferase             | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.05 | 0.10 | 0.15 | 0.14 |
| 13         | 890,427  | G>C      | fhuC                       | G217A(GGA>GCA)      | Iron(3+)hydroxamate import ATPbinding protein FhuC                 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.24 | 0.12 | 0.03 |
| 14         | 955,123  | C>A      | PROKKA_00981<              | G119C(GGC>TGC)      | putative HTHtype transcriptional regulator                         | 0.00 | 0.00 | 0.05 | 0.20 | 0.26 | 0.10 | 0.26 | 0.60 | 0.66 | 0.88 | 0.97 |

Note that the frequency of each mutation can be reported as either a percentage out of 100 (ex. 79.42%) or as a number betweeThe two-step method was originally devised by Katya Kosheleva in 2012 to model the evolution of yeast populations, and has since been modified to accomodate a wider array of experimental designs. Since each frequency measurement is analagous to the probability of a mutation being detected at a given timepoint, we can use the binomial distribution to test whether two sets of frequency measurements represent the same underlying series. There are a few key assumptions that must be made:n 0 and 1 (ex. 0.7942). Regardless of which format is contained in the input table, the scripts will convert the frequency values so they fit in the range $[0,1]$.
