#!/usr/bin/env python3
'''
Map conservation data to an alignment
'''

import sys

import bioinfo


def map_cons(headers, seqs, data):
    '''
    Map a conserved position to the alignment
    '''
    prot = [header for header in headers if data[0] in header][0]
    pos_rat = data[1:]
    ali_data = [0.0]*len(seqs[0])
    for posi in range(len(pos_rat[0])):
        pos = map_pos(pos_rat[0][posi], seqs[headers.index(prot)])
        ali_data[pos] = pos_rat[1][posi]
    return ali_data


def map_pos(pos, sequence):
    '''
    Map a position in a protein to the position in the alignment
    input: pos1 = 1
    output: pos1 = 0
    '''
    p = 0
    for i in range(len(sequence)):
        if sequence[i] != '-':
            p += 1
        if p == pos:
            return i
    err = 'E: position {} not found in {}\n'.format(pos, sequence[:20])
    sys.stderr.write(err)
    return False


def main(fasta_name, data_files):
    '''
    Map conservation data to an alignment
    '''
    heads, seqs = bioinfo.read_fasta(fasta_name)
    consnames = list()
    ali_cons = list()
    for d_file in data_files:
        cons = read_data(d_file)
        if cons is False:
            return False
        ali_cons.append(map_cons(heads, seqs, cons))
        consnames.append(cons[0])
    for i in range(len(seqs)):
        print('{}\t{}'.format('ali:' + heads[i],
                              '\t'.join([c for c in seqs[i]])))
    for i in range(len(ali_cons)):
        print('{}\t{}'.format('cons:' + consnames[i],
                              '\t'.join([str(r) for r in ali_cons[i]])))


def read_data(filename):
    '''
    Read conservation data
    '''
    with open(filename) as infile:
        raw = infile.read()
    lines = raw.split('\n')
    if lines[0][0] != '#':
        sys.stderr.write('E: the file {} has no reference\n'.format(filename))
        return False
    protname = lines[0][1:].strip()
    lines = [l for l in lines[1:] if len(l) > 0]
    data = [line.split('\t') for line in lines]
    positions = [int(dat[0]) for dat in data]
    ratios = [float(dat[3]) for dat in data]
    return (protname, positions, ratios)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        USAGE = 'Usage: {} <fasta alignment> <data files>\n'.format(sys.argv[0])
        sys.stderr.write(USAGE)
        sys.exit(1)
    if main(sys.argv[1], sys.argv[2:]) is False:
        sys.exit(1)
