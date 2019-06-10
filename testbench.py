import os

import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-p', '--path', required=True, type=str, default='', help="Path to a folder with Multiple Alignment File (MAF)")
parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File prefix")
args = parser.parse_args()

os.system("python ./main.py -m " + args.path + " -o " + args.output)
print("LSE of non corrected matrix")
os.system("python ./LSEdiag.py -d " + args.output + "_no_corr.csv" + "| cut -f2 | awk '{n += $1}; END{print n}'")
print("LSE of corrected matrix")
os.system("python ./LSEdiag.py -d " + args.output + "_corr.csv" + "| cut -f2 | awk '{n += $1}; END{print n}'")