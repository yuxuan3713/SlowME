import os
import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-m', '--maf', required=True, type=str, default='', help="Path to a Multiple Alignment File (MAF)")
parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File prefix")
# parser.add_argument('-i', '--id', required=False, type=str, default='', help="ID files")
args = parser.parse_args()

all_not_no_corr = args.output + "_all_not" + "_no_corr.csv"
all_not_corr = args.output + "_all_not" + "_corr.csv"
pair_not_no_corr = args.output + "_pair_not" + "_no_corr.csv"
pair_not_corr = args.output + "_pair_not" + "_corr.csv"
mismatch_no_corr = args.output + "_mismatch" + "_no_corr.csv"
mismatch_corr = args.output + "_mismatch" + "_corr.csv"


os.system("python ./main.py -m %s -o %s " % (args.maf, args.output))


print("LSE of non corrected matrix [pair not non-corr]")
os.system("python ./LSEdiag.py -d %s | cut -f2 | awk '{n += $1}; END{print n}'" % pair_not_no_corr)
print("LSE of non corrected matrix [pair not corr]")
os.system("python ./LSEdiag.py -d %s | cut -f2 | awk '{n += $1}; END{print n}'" % pair_not_corr)
print("LSE of non corrected matrix [all not non-corr]")
os.system("python ./LSEdiag.py -d %s | cut -f2 | awk '{n += $1}; END{print n}'" % all_not_no_corr)
print("LSE of non corrected matrix [all not corr]")
os.system("python ./LSEdiag.py -d %s | cut -f2 | awk '{n += $1}; END{print n}'" % all_not_corr)
print("LSE of non corrected matrix [mismatch non-corr]")
os.system("python ./LSEdiag.py -d %s | cut -f2 | awk '{n += $1}; END{print n}'" % mismatch_no_corr)
print("LSE of non corrected matrix [mismatch corr]")
os.system("python ./LSEdiag.py -d %s | cut -f2 | awk '{n += $1}; END{print n}'" % mismatch_corr)