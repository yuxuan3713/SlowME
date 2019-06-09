from Bio import AlignIO
import numpy as np
import scipy as sp
from itertools import combinations
from tool_func import mutation_count, segment_mut_count, get_gamma_param, to_phylip_format, to_csv_format
from math import log
import sys

import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-m', '--maf', required=True, type=str, default='', help="Multiple Alignment File (MAF)")
parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Distance Matrix")
args = parser.parse_args()

# if args.maf == '':
#     filePath = '/Users/liuyuxuan/Downloads/208dataset/dog_chrX.maf'
# else:
#     filePath = args.maf

filePath = args.maf
SEQ_LIMIT = 1000


mat = dict()
allID = set()
alignment_count = 0
seq_count = 0


# one pass to get all id
print("Scanning for all IDs...")
progress_count = 0
for alignment_block in AlignIO.parse(filePath, 'maf'):
    alignment_count += 1
    for seqrec in alignment_block:
        allID.add(seqrec.id.split('.')[0])  # only consider species
        seq_count += 1

# initialize a matrix
allID = sorted(allID)
for id1 in allID:
    mat[id1] = dict()
    for id2 in allID:
        mat[id1][id2] = dict()
        mat[id1][id2]['total_mut'] = 0
        mat[id1][id2]['total_len'] = 0
        mat[id1][id2]['D'] = 0
        mat[id1][id2]['total_seq'] = 0
        mat[id1][id2]['data_point'] = list()
        mat[id1][id2]['gamma'] = dict()
        mat[id1][id2]['prob_mut'] = 0

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
        print('Progress [%d%%]' % progress)
        sys.stdout.flush()

    comb = combinations(list(alignment_block), 2)
    valid_block += 1
    valid_seq += alignment_block.get_alignment_length()
    for seqrec1, seqrec2 in comb:
        if seqrec1.annotations['size'] < SEQ_LIMIT or seqrec2.annotations['size'] < SEQ_LIMIT:
            continue
        # start processing on the alignment
        valid_block += 1
        valid_seq += alignment_block.get_alignment_length()
        pair_id = sorted([seqrec1.id.split('.')[0], seqrec2.id.split('.')[0]])
        id1 = pair_id[0]
        id2 = pair_id[1]
        mat[id1][id2]['total_seq'] += 1
        mut_count, seq_len = mutation_count(seqrec1.seq, seqrec2.seq)
        mat[id1][id2]['total_mut'] += mut_count
        mat[id1][id2]['total_len'] += seq_len
        mat[id1][id2]['data_point'].extend(segment_mut_count(seqrec1.seq, seqrec2.seq, SEQ_LIMIT))


# process data in the matrix
for i in range(len(allID)):
    for j in range(i+1, len(allID)):
        id1 = allID[i]
        id2 = allID[j]
        if mat[id1][id2]['total_len'] == 0:
            allID.remove(id1)
            allID.remove(id2)
            del mat[id1][id2]
            del mat[id2][id1]
        else:
            mat[id1][id2]['prob_mut'] = mat[id1][id2]['total_mut'] / mat[id1][id2]['total_len']
            mat[id1][id2]['d'] = -3/4 * log(1 - 4/3 * mat[id1][id2]['prob_mut'])
            mat[id1][id2]['data_point'] = np.array([x / SEQ_LIMIT for x in mat[id1][id2]['data_point']])
            mat[id1][id2]['gamma']['param'] = get_gamma_param(mat[id1][id2]['data_point'])
            # mat[id1][id2]['gamma']['alpha'] = mat[id1][id2]['gamma']['param'][0]
            mat[id1][id2]['gamma']['alpha'] = 1 / np.var(mat[id1][id2]['data_point'] / np.mean(mat[id1][id2]['data_point']))

# no correction: from D to matrix
mat_no_corr = list()
for id1 in allID:
    temp = list()
    for id2 in allID:
        temp.append(mat[id1][id2]['prob_mut'])
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

print("end of program\n")

for id1 in allID:
    for id2 in allID:
        if 'alpha' in mat[id1][id2]['gamma']:
            alpha = mat[id1][id2]['gamma']['alpha']
            d = mat[id1][id2]['d']
            t = 3 / 4 * alpha * ((1 - 4 / 3 * d) ** (-1 / alpha) - 1)
            print(mat[id1][id2]['gamma']['param'],t)


from sys import stdout as outfile
to_phylip_format(allID, mat_no_corr, outfile)
to_phylip_format(allID, mat_corr, outfile)

if args.output == 'stdout':
    to_csv_format(allID, mat_no_corr, outfile)
    to_csv_format(allID, mat_corr, outfile)
else:
    print("Printing to file...")
    outfile1 = open(args.output + "_no_corr.phylip", 'w')
    to_phylip_format(allID, mat_no_corr, outfile1)
    outfile1.close()
    outfile2 = open(args.output + "_corr.phylip", 'w')
    to_phylip_format(allID, mat_corr, outfile2)
    outfile2.close()
    np.save(args.output + "_no_corr.npy", mat_no_corr)
    np.save(args.output + "_corr.npy", mat_corr)

print("All id(s):")
print(allID)