# ECE-208-project

[main.py](main.py) will take as input the following items:
* `-m`: The path to a [maf file](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) containing multiple sequence alignment
* `-o`: [optional] Output file path and name. Name is used as prefix of filename for all outfiles. Required for output file.
* `-l`: [optional] The window size used for gamma distribution model
* `-v`: [optional] Verbose mode flag. Pass False to suppress explicit result output. Default as True.
* `-i`: [optional] The path to a ID file that corresponding to the MAF input file.
* `-s`: [optional] Mismatch strategy selection. Default as try all strategy.

[testBench.py](testBench.py) will take as input the following items:
* `-m`: The path to a [maf file](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) containing multiple sequence alignment
* `-o`: [optional] Output file path and name. Name is used as prefix of filename for all outfiles. Required for output file.

