# ECE-208-project: SlowME
This repository is a tool that gives estimate of phylogenetic tree based on multiple alignment data. The name SlowME is a joke that we found FastME to be a powerful and great tool in the background of computational biology and that our program takes about 40 hours running 8GB maf file on a single CPU. The program scans through the maf file, retrives all species IDs (if you don't have it handy), takes long segment of aligned sequence and estimate a gamma distribution of the local mutation rate. It computes a distance matrix of all pair of species in the maf file. FastME takes distance matrix and gives phylogenetic tree.


### Acknowlegement:
* Professor Siavash Mirarab's help on theory and direction of the project.
* LSE script from Metin Balaban, [LSEdiag](https://github.com/balabanmetin/LSEdiag)
* The great tool [FastME](http://www.atgc-montpellier.fr/fastme/) that inspires us on both our project content and our project name.
* Multiple Alignment dataset from [UCSC Genome Browser](https://genome.ucsc.edu)

### Dependency
* [FastME](http://www.atgc-montpellier.fr/fastme/)
* [TreeSwift](https://github.com/niemasd/TreeSwift)
* [Biopython](https://biopython.org)
* [LSEdiag](https://github.com/balabanmetin/LSEdiag)

### Usage
[main.py](main.py) is the python script that computes distance matrix. It will take as input the following items:
* `-m`: The path to a [maf file](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) containing multiple sequence alignment
* `-o`: [optional] Output file path and name. Name is used as prefix of filename for all outfiles. Required for output file.
* `-l`: [optional] The window size used for gamma distribution model
* `-v`: [optional] Verbose mode flag. Pass False to suppress explicit result output. Default as True.
* `-i`: [optional] The path to a ID file that corresponding to the MAF input file.
* `-s`: [optional] Mismatch strategy selection. Default as try all strategy.

[testbench.py](testbench.py) combines calling main.py and FastME. It will take as input the following items:
* `-m`: The path to a [maf file](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) containing multiple sequence alignment
* `-o`: [optional] Output file path and name. Name is used as prefix of filename for all outfiles. Required for output file.

[LSEdiag.py](LSEdiag.py) is the slightly modified version from [LSEdiag](https://github.com/balabanmetin/LSEdiag)

* Script for getting LSE:
```
python ~/LSEdiag/LSEdiag.py -d test.mat -t  test.nwk | cut -f2 | awk '{n += $1}; END{print n}'
```

[test.sh](test.sh) is a simple shell script that helps run LSEdiag on all csv results.

### Included Sample Result
* csv file of distance matrix on [Dog, Mouse, Rat, Human 4-way alignment](http://hgdownload.soe.ucsc.edu/goldenPath/canFam2/multiz4way/)
* csv file of distance matrix on select [Human 20 way multiple alignment](http://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz20way/README.txt)
