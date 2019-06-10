from Bio import AlignIO
import numpy as np
# import scipy as sp
# import subprocess
from itertools import combinations
from tool_func import get_gamma_param, to_phylip_format, to_csv_format
from tool_func import mutation_count_pair_not, segment_mut_count_pair_not
from tool_func import mutation_count_all_not, segment_mut_count_all_not, process_block_all_not
from tool_func import mutation_count_mismatch, segment_mut_count_mismatch
from math import log
import sys

import argparse

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-m', '--maf', required=True, type=str, help="Multiple Alignment File (MAF)")
parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output Distance Matrix")
parser.add_argument('-l', '--limit', required=False, type=int, default=1000, help="Sequence length limit")
parser.add_argument('-s', '--strategy', required=False, type=int, default=0, help="Indel Strategy")
parser.add_argument('-i', '--id', required=False, type=str, default='', help="ID files")
parser.add_argument('-v', '--verbose', required=False, type=bool, default=True, help="Verbose mode")
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

        mat[id1][id2]['pair_not'] = dict()
        mat[id1][id2]['pair_not']['mutation'] = 0                   # count of mutation
        mat[id1][id2]['pair_not']['total'] = 0                      # total length of scan sequence
        mat[id1][id2]['pair_not']['d'] = 0                          # humming distance
        mat[id1][id2]['pair_not']['t_non_corr'] = 0                 # non corrected phylogenetic dist
        mat[id1][id2]['pair_not']['t_corr'] = 0                     # corrected phylogenetic distance
        mat[id1][id2]['pair_not']['gamma'] = dict()
        mat[id1][id2]['pair_not']['gamma']['data_point'] = list()   # gamma data point

        mat[id1][id2]['all_not'] = dict()
        mat[id1][id2]['all_not']['mutation'] = 0                    # count of mutation
        mat[id1][id2]['all_not']['total'] = 0                       # total length of scan sequence
        mat[id1][id2]['all_not']['d'] = 0                           # humming distance
        mat[id1][id2]['all_not']['t_non_corr'] = 0                  # non corrected phylogenetic dist
        mat[id1][id2]['all_not']['t_corr'] = 0                      # corrected phylogenetic distance
        mat[id1][id2]['all_not']['gamma'] = dict()
        mat[id1][id2]['all_not']['gamma']['data_point'] = list()    # gamma data point

        mat[id1][id2]['mismatch'] = dict()
        mat[id1][id2]['mismatch']['mutation'] = 0                   # count of mutation
        mat[id1][id2]['mismatch']['total'] = 0                      # total length of scan sequence
        mat[id1][id2]['mismatch']['d'] = 0                          # humming distance
        mat[id1][id2]['mismatch']['t_non_corr'] = 0                 # non corrected phylogenetic dist
        mat[id1][id2]['mismatch']['t_corr'] = 0                     # corrected phylogenetic distance
        mat[id1][id2]['mismatch']['gamma'] = dict()
        mat[id1][id2]['mismatch']['gamma']['data_point'] = list()   # gamma data point

        # mat[id1][id2]['total_mut'] = 0
        # mat[id1][id2]['total_len'] = 0
        # mat[id1][id2]['D'] = 0
        # mat[id1][id2]['total_seq'] = 0
        # mat[id1][id2]['data_point'] = list()
        # mat[id1][id2]['gamma'] = dict()
        # mat[id1][id2]['humming'] = 0

# process raw data
print("Processing raw data...")
progress_count = 0
progress = 0
valid_block = 0
valid_seq = 0
for alignment_block in AlignIO.parse(filePath, 'maf'):

    # update progress counter
    progress_count += 1
    if round(progress_count / alignment_count * 100) > progress:
        progress = round(progress_count / alignment_count * 100)
        # print('Progress [%d%%]' % progress)
        print('Progress [%d%%]\r' % progress, end="")
        sys.stdout.flush()

    # statistics count valid block
    # valid_block += 1
    # valid_seq += len(alignment_block)

    # strategy: all not
    # all_seqrec = [seqrec for seqrec in alignment_block]
    # TODO: process all sequence, remove indel column
    all_seqrec = process_block_all_not(alignment_block)
    comb = combinations(all_seqrec, 2)
    for seqrec1, seqrec2 in comb:
        pair_id = sorted([seqrec1.id.split('.')[0], seqrec2.id.split('.')[0]])
        id1 = pair_id[0]
        id2 = pair_id[1]

        mut_count, seq_len = mutation_count_mismatch(seqrec1.seq, seqrec2.seq)

        mat[id1][id2]['all_not']['mutation'] += mut_count
        mat[id1][id2]['all_not']['total'] += seq_len

        if not (len(seqrec1) < SEQ_LIMIT or len(seqrec2) < SEQ_LIMIT):
            mat[id1][id2]['all_not']['gamma']['data_point'].extend(
                segment_mut_count_mismatch(seqrec1.seq, seqrec2.seq, SEQ_LIMIT))

    # strategy: mismatch and pair not
    comb = combinations(list(alignment_block), 2)
    for seqrec1, seqrec2 in comb:

        # start processing on the alignment
        valid_block += 1
        valid_seq += alignment_block.get_alignment_length()
        pair_id = sorted([seqrec1.id.split('.')[0], seqrec2.id.split('.')[0]])
        id1 = pair_id[0]
        id2 = pair_id[1]
        # mat[id1][id2]['total_seq'] += 1  # count total sequence pairs

        # indel strategy 1 : count as not for pair
        mut_count, seq_len = mutation_count_pair_not(seqrec1.seq, seqrec2.seq)

        mat[id1][id2]['pair_not']['mutation'] += mut_count
        mat[id1][id2]['pair_not']['total'] += seq_len

        if not (seqrec1.annotations['size'] < SEQ_LIMIT or seqrec2.annotations['size'] < SEQ_LIMIT):
            mat[id1][id2]['pair_not']['gamma']['data_point'].extend(segment_mut_count_pair_not(seqrec1.seq, seqrec2.seq, SEQ_LIMIT))

        # indel strategy: count as a mismatch
        mut_count, seq_len = mutation_count_mismatch(seqrec1.seq, seqrec2.seq)

        mat[id1][id2]['mismatch']['mutation'] += mut_count
        mat[id1][id2]['mismatch']['total'] += seq_len

        if not (seqrec1.annotations['size'] < SEQ_LIMIT or seqrec2.annotations['size'] < SEQ_LIMIT):
            mat[id1][id2]['mismatch']['gamma']['data_point'].extend(segment_mut_count_mismatch(seqrec1.seq, seqrec2.seq, SEQ_LIMIT))


# process data in the matrix
invalid_list = list()
strategy_list = ['pair_not', 'all_not', 'mismatch']
for i in range(len(allID)):
    for j in range(i+1, len(allID)):
        for strategy in strategy_list:
            id1 = allID[i]
            id2 = allID[j]
            if mat[id1][id2][strategy]['total'] == 0:
                invalid_list.append(id1)
                invalid_list.append(id2)
                print("***** Warning: distance matrix not full *****")
            else:
                # non corrected phylogenetic distance
                mat[id1][id2][strategy]['d'] = mat[id1][id2][strategy]['mutation'] / mat[id1][id2][strategy]['total']
                mat[id1][id2][strategy]['t_non_corr'] = -3/4 * log(1 - 4/3 * mat[id1][id2][strategy]['d'])

                if len(mat[id1][id2][strategy]['gamma']['data_point']) > 0:
                    mat[id1][id2][strategy]['gamma']['data_point'] = np.array([x / SEQ_LIMIT for x in mat[id1][id2][strategy]['gamma']['data_point']])
                    mat[id1][id2][strategy]['gamma']['alpha'] = 1 / np.var(mat[id1][id2][strategy]['gamma']['data_point'] / np.mean(mat[id1][id2][strategy]['gamma']['data_point']))
                    alpha = mat[id1][id2][strategy]['gamma']['alpha']
                    d = mat[id1][id2][strategy]['d']
                    mat[id1][id2][strategy]['t_corr'] = 3 / 4 * alpha * ((1 - 4 / 3 * d) ** (-1 / alpha) - 1)
                else:
                    # no correction data point, no seq rec above SEQ_LIMIT
                    print("***** Warning: No correction data point for %s and %s *****", id1, id2)

                # mat[id1][id2][strategy]['gamma']['param'] = get_gamma_param(mat[id1][id2][strategy]['data_point'])
                # mat[id1][id2]['gamma']['alpha'] = mat[id1][id2]['gamma']['param'][0] # deprecated gamma fitting

# # remove invalid ID entries
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

mat_non_corr = dict()
mat_corr = dict()
for strategy in strategy_list:
    # no correction distance matrix
    mat_non_corr[strategy] = list()
    for id1 in allID:
        temp = list()
        for id2 in allID:
            temp.append(mat[id1][id2][strategy]['d'])
        mat_non_corr[strategy].append(temp)
    for i in range(len(allID)):
        mat_non_corr[strategy][i][i] = 0
    mat_non_corr[strategy] = np.array(mat_non_corr[strategy]) + np.array(mat_non_corr[strategy]).transpose()

    # with correction: from alpha to distance matrix
    mat_corr[strategy] = list()
    for id1 in allID:
        temp = list()
        for id2 in allID:
            temp.append(mat[id1][id2][strategy]['t_corr'])
        mat_corr[strategy].append(temp)
    for i in range(len(allID)):
        mat_corr[strategy][i][i] = 0
    mat_corr[strategy] = np.array(mat_corr[strategy]) + np.array(mat_corr[strategy]).transpose()


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

if args.verbose:
    from sys import stdout as outfile
    print("********************Printing result********************\n")
    print("[All id(s)]")
    print(allID)
    for strategy in strategy_list:
        print("[Non-corrected matrix: " + strategy + "]")
        to_phylip_format(allID, mat_non_corr[strategy], outfile)
        print("\n[Corrected matrix: " + strategy + "]")
        to_phylip_format(allID, mat_corr[strategy], outfile)
    print("\n*********************End of Result*********************")

if args.output != 'stdout':
    print("Printing to file...")
    # to_phylip_format(allID, mat_no_corr, outfile)
    # to_phylip_format(allID, mat_corr, outfile)

    all_not_no_corr = args.output + "_all_not" + "_no_corr.csv"
    outfile = open(all_not_no_corr, 'w')
    to_csv_format(allID, mat_non_corr['all_not'], outfile)
    outfile.close()

    all_not_corr = args.output + "_all_not" + "_corr.csv"
    outfile = open(all_not_corr, 'w')
    to_csv_format(allID, mat_corr['all_not'], outfile)
    outfile.close()

    pair_not_no_corr = args.output + "_pair_not" + "_no_corr.csv"
    outfile = open(pair_not_no_corr, 'w')
    to_csv_format(allID, mat_non_corr['pair_not'], outfile)
    outfile.close()

    pair_not_corr = args.output + "_pair_not" + "_corr.csv"
    outfile = open(pair_not_corr, 'w')
    to_csv_format(allID, mat_corr['pair_not'], outfile)
    outfile.close()

    mismatch_no_corr = args.output + "_mismatch" + "_no_corr.csv"
    outfile = open(mismatch_no_corr, 'w')
    to_csv_format(allID, mat_non_corr['mismatch'], outfile)
    outfile.close()

    mismatch_corr = args.output + "_mismatch" + "_corr.csv"
    outfile = open(mismatch_corr, 'w')
    to_csv_format(allID, mat_corr['mismatch'], outfile)
    outfile.close()

    # np.save(args.output + "_no_corr.npy", mat_no_corr)
    # np.save(args.output + "_corr.npy", mat_corr)

# s = ["fastme", "-i", dist_phy.name, "-o", tree_fp]
# subprocess.call(s, stdout = nldef, stderr = nldef)

print("[Program Finished]\n")
