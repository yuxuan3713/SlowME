from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
aln = AlignIO.parse("/Users/liuyuxuan/Downloads/208dataset/dog_chrX.maf", "maf")
AlignIO.write(aln, "/Users/liuyuxuan/Downloads/208dataset/dog_chrX.fa", "fasta")
# calculator = DistanceCalculator('identity')
# dm = calculator.get_distance(aln)
# print