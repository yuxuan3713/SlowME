import scipy as sp
import scipy.stats


# return a tuple (mutation count, length)
# strategy: count as not for pair
def mutation_count_pair_not(r, s):
    assert len(r) == len(s), "sequence not aligned"
    count = 0
    length = len(r)
    for i in range(len(r)):
        if r[i] == '-' or s[i] == '-':
            length -= 1
            continue
        elif r[i].upper() != s[i].upper():
            count += 1
    return count, length


# return a tuple (mutation count, length)
# strategy: count as mismatch
def mutation_count_mismatch(r, s):
    assert len(r) == len(s), "sequence not aligned"
    count = 0
    length = len(r)
    for i in range(len(r)):
        if r[i].upper() != s[i].upper():
            count += 1
    return count, length


def segment_mut_count_pair_not(r, s, l):
    assert len(r) == len(s)
    assert len(r) > l
    result = list()
    for i in range(int(len(r) / l)):
        first_seq = r[l * i: l * (i + 1) - 1]
        second_seq = s[l * i: l * (i + 1) - 1]
        result.append(mutation_count_pair_not(first_seq, second_seq)[0])
    return result


def segment_mut_count_mismatch(r, s, l):
    assert len(r) == len(s)
    assert len(r) > l
    result = list()
    for i in range(int(len(r) / l)):
        first_seq = r[l * i: l * (i + 1) - 1]
        second_seq = s[l * i: l * (i + 1) - 1]
        result.append(mutation_count_mismatch(first_seq, second_seq)[0])
    return result


# TODO: all not
def mutation_count_all_not(block, s):
    for seqrec in block:
        print()
    return 0


# TODO: all not
def segment_mut_count_all_not(r, s, l):
    return 0


# take a list of data point and return the variance alpha of a gamma distribution
# option to return mse of the fit?
def get_gamma_param(data):
    data = data / sp.mean(data)
    dist = getattr(scipy.stats, 'gamma')
    param = dist.fit(data)
    return param


def to_phylip_format(ids, data, output):
    id_len = max([len(i) for i in ids])
    # data_len = 8
    output.write('%d\n' % len(ids))
    fstring = '  '.join(['%.8f' for i in range(len(ids))]) + '\n'
    for i in range(len(ids)):
        output.write('%-*s' % (id_len+2, ids[i]))
        output.write(fstring % tuple(dist for dist in data[i]))


def to_csv_format(ids, data, output):
    fstring = ''.join(['\t%s' for i in ids]) + '\n'
    output.write(fstring % tuple(i for i in ids))
    fstring = ''.join(['%.8f\t' for i in ids]) + '\n'
    for i in range(len(ids)):
        output.write('%s\t' % ids[i])
        output.write(fstring % tuple(dist for dist in data[i]))
