#!/usr/bin/env python3
'''
Map conservation data to an alignment
'''

import sys

import bioinfo


def map_cons(headers, seqs, ali_data, data):
    '''
    Map a conserved position to the alignment
    '''
    prot = [header for header in headers if data[0] in head]
    for dat in range(len(data)):
        pos = map_posdata(data[dat][0]-1, seqs[headers.index(prot)])
        ali_data[pos] = data[dat][3]
    return ali_data

def map_pos(pos, protein):
    '''
    Map a position in a protein to the position in the alignment
    Pos 1 = 0
    '''
    pass


def main(fasta_name, data_files):
    '''
    Map conservation data to an alignment
    '''
    heads, seqs = bioinfo.read_fasta(fasta_name)
    data = list()
    ali_cons = list()
    for d_file in data_files:
        cons = read_data(d_file)
        if cons == -1:
            return -1
        data.append(read_data(cons))
        ali_cols.append([0.0]*len(seqs[0]))
    
    
    # put data in correct positions
    # print alignment with data
    # option to only show positions with any data > 


def read_data(filename):
    '''
    Read conservation data
    '''
    with open(filename) as infile:
        raw = infile.read()
    lines = raw.split('\n')
    if lines[0][0] != '#':
        sys.stderr.write('E: the file {} has no reference\n'.format(filename))
        return -1
    protname = lines[0][1:].strip()
    lines = [l for l in lines[1:] if len(l) > 0]
    data = [line.split('\t') for line in lines]
    positions = [int(dat[0]) for dat in data]
    ratios = [float(dat[3]) for dat in data]
    return (protname, positions, ratios)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.stderr.write('Usage: {} <fasta alignment> <data files>\n'.format(sys.argv[0]))
        sys.exit(1)
    if main(sys.argv[1], sys.argv[2:]) == -1:
        sys.exit(1)
    
