from Bio import AlignIO
import numpy as np
# import scipy as sp
# import subprocess
from itertools import combinations
from tool_func import mutation_count_pair_not, segment_mut_count_pair_not, get_gamma_param, to_phylip_format, to_csv_format
from math import log
import sys

import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-m', '--maf', required=True, type=str, help="Multiple Alignment File (MAF)")
parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Distance Matrix")
parser.add_argument('-l', '--limit', required=False, type=int, default=1000, help="Sequence length limit")
parser.add_argument('-s', '--strategy', required=False, type=int, default=0, help="Indel Strategy")
parser.add_argument('-i', '--id', required=False, type=str, default='', help="ID files")
args = parser.parse_args()

filePath = args.maf
SEQ_LIMIT = args.limit


mat = dict()
allID = set()
alignment_count = 0
seq_count = 0


# one pass to get all id
if len(args.id) == 0:
    print("Scanning for all IDs...")
    progress_count = 0
    for alignment_block in AlignIO.parse(filePath, 'maf'):
        alignment_count += 1
        for seqrec in alignment_block:
            allID.add(seqrec.id.split('.')[0])  # only consider species
            seq_count += 1
else:
    print("Reading from ID file...")
    allID = sorted([line.rstrip('\n') for line in open(args.id)])

print("***** All ID *****")
for i in allID:
    print(i)
print("*** End of IDs ***")

# initialize a matrix
allID = sorted(allID)
for id1 in allID:
    mat[id1] = dict()
    for id2 in allID:
        mat[id1][id2] = dict()
        mat[id1][id2]['total_mut_pair_not'] = 0
        mat[id1][id2]['total_len'] = 0
        mat[id1][id2]['D'] = 0
        mat[id1][id2]['total_seq'] = 0
        mat[id1][id2]['data_point'] = list()
        mat[id1][id2]['gamma'] = dict()
        mat[id1][id2]['humming'] = 0

# process raw data
print("Processing raw data...")
progress_count = 0
progress = 0
valid_block = 0
valid_seq = 0
for alignment_block in AlignIO.parse(filePath, 'maf'):
    progress_count += 1
    if round(progress_count / alignment_count * 100) > progress:
        progress = round(progress_count / alignment_count * 100)
        # print('Progress [%d%%]' % progress)
        print('Progress [%d%%]\r' % progress, end="")
        sys.stdout.flush()

    comb = combinations(list(alignment_block), 2)
    # valid_block += 1
    # valid_seq += len(alignment_block)
    for seqrec1, seqrec2 in comb:
        if seqrec1.annotations['size'] < SEQ_LIMIT or seqrec2.annotations['size'] < SEQ_LIMIT:
            continue
        # start processing on the alignment
        valid_block += 1
        valid_seq += alignment_block.get_alignment_length()
        pair_id = sorted([seqrec1.id.split('.')[0], seqrec2.id.split('.')[0]])
        id1 = pair_id[0]
        id2 = pair_id[1]
        mat[id1][id2]['total_seq'] += 1  # count total sequence pairs

        # indel strategy 1 : count as not for pair
        mut_count, seq_len = mutation_count_pair_not(seqrec1.seq, seqrec2.seq)
        mat[id1][id2]['total_mut'] += mut_count
        mat[id1][id2]['total_len'] += seq_len
        mat[id1][id2]['data_point'].extend(segment_mut_count_pair_not(seqrec1.seq, seqrec2.seq, SEQ_LIMIT))
        # indel strategy 2 : count as not for all


        # indel strategy 3 : count as a mismatch




# process data in the matrix
invalid_list = list()
for i in range(len(allID)):
    for j in range(i+1, len(allID)):
        id1 = allID[i]
        id2 = allID[j]
        if mat[id1][id2]['total_len'] == 0:
            invalid_list.append(id1)
            invalid_list.append(id2)
            print("***** Warning: distance matrix not full *****")
        else:
            mat[id1][id2]['humming'] = mat[id1][id2]['total_mut'] / mat[id1][id2]['total_len']
            mat[id1][id2]['d'] = -3/4 * log(1 - 4/3 * mat[id1][id2]['humming'])
            mat[id1][id2]['data_point'] = np.array([x / SEQ_LIMIT for x in mat[id1][id2]['data_point']])
            mat[id1][id2]['gamma']['param'] = get_gamma_param(mat[id1][id2]['data_point'])
            # mat[id1][id2]['gamma']['alpha'] = mat[id1][id2]['gamma']['param'][0]
            mat[id1][id2]['gamma']['alpha'] = 1 / np.var(mat[id1][id2]['data_point'] / np.mean(mat[id1][id2]['data_point']))

# for id1 in invalid_list:
#     for id2 in invalid_list:
#         if id1 in allID:
#             allID.remove(id1)
#             del mat[id1]
#         if id2 in allID:
#             allID.remove(id2)
#             del mat[id2]
if len(invalid_list) != 0:
    print("Empty entry with following IDs: ", invalid_list)

# no correction: from D to matrix
mat_no_corr = list()
for id1 in allID:
    temp = list()
    for id2 in allID:
        temp.append(mat[id1][id2]['humming'])
    mat_no_corr.append(temp)
for i in range(len(allID)):
    mat_no_corr[i][i] = 0
mat_no_corr = np.array(mat_no_corr) + np.array(mat_no_corr).transpose()

# with correction: from alpha to matrix
mat_corr = list()
for id1 in allID:
    temp = list()
    for id2 in allID:
        if 'alpha' in mat[id1][id2]['gamma']:
            alpha = mat[id1][id2]['gamma']['alpha']
            d = mat[id1][id2]['d']
            t = 3 / 4 * alpha * ((1 - 4 / 3 * d) ** (-1/alpha) - 1)
            temp.append(t)
        else:
            temp.append(0)
    mat_corr.append(temp)
for i in range(len(allID)):
    mat_corr[i][i] = 0
mat_corr = np.array(mat_corr) + np.array(mat_corr).transpose()


# print result
# print('alignment count:', alignment_count)
# print('seq count:', seq_count)
#
# print('valid alignment count pair:', valid_block)
# print('valid seq count pair:', valid_seq)
#
# print('Printing non-corrected matrix:')
# print(sorted(allID))
# print(mat_no_corr)
#
# print('Printing corrected matrix:')
# print(sorted(allID))
# print(mat_corr)

print("[Processing Complete]\n")

# for id1 in allID:
#     for id2 in allID:
#         if 'alpha' in mat[id1][id2]['gamma']:
#             alpha = mat[id1][id2]['gamma']['alpha']
#             d = mat[id1][id2]['d']
#             t = 3 / 4 * alpha * ((1 - 4 / 3 * d) ** (-1 / alpha) - 1)
#             print(mat[id1][id2]['gamma']['param'],t)


from sys import stdout as outfile
print("********************Printing result********************")
print("[Non-corrected matrix]")
to_phylip_format(allID, mat_no_corr, outfile)
print("\n[Corrected matrix]")
to_phylip_format(allID, mat_corr, outfile)

print("\n[All id(s)]")
print(allID)
print("\n*********************End of Result*********************")

if args.output != 'stdout':

    # to_phylip_format(allID, mat_no_corr, outfile)
    # to_phylip_format(allID, mat_corr, outfile)

    print("Printing to file...")
    outfile1 = open(args.output + "_no_corr.csv", 'w')
    to_csv_format(allID, mat_no_corr, outfile1)
    outfile1.close()
    outfile2 = open(args.output + "_corr.csv", 'w')
    to_csv_format(allID, mat_corr, outfile2)
    outfile2.close()
    # np.save(args.output + "_no_corr.npy", mat_no_corr)
    # np.save(args.output + "_corr.npy", mat_corr)

# s = ["fastme", "-i", dist_phy.name, "-o", tree_fp]
# subprocess.call(s, stdout = nldef, stderr = nldef)


