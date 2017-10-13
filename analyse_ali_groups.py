#!/usr/bin/env python3
'''
Analyse a protein sequence alignment based on properties rather
than exact residues.
'''

import sys

import bioinfo

GROUPS = {'I':'B'}

def conservation(sequences):
    '''
    Analyse an alignment based on conservation.
    '''
    best_cons = [0] * len(sequences[0])
    nseqs = len(sequences)
    for i in range(len(sequences[0])):
        pos_res = tuple(seq[i] for seq in sequences)
        res_counts = tuple((res, pos_res.count(res)/nseqs) for res in set(pos_res))
        counts = tuple(count[1] for count in res_counts)
        best_cons[i] = res_counts[counts.index(max(counts))]
    return best_cons


def main(filename):
    '''
    Read a FASTA alignment and perform a conservation analysis
    using groups of residues
    '''
    sequences = bioinfo.read_fasta(filename)[1]

    sequences = transform_prop(sequences)
    cons = conservation(sequences)
    for con in cons:
        print('{}\t{:.2f}'.format(con[0], con[1]))


def transform_prop(sequences):
    '''
    Input: sequences - list of protein sequences
    Return: the same alignment with the relevant residues replaced
    by their group names
    Transform an alignment by grouping residues:
    I - ILV
    D - DE
    K - KR
    '''
    for i in range(len(sequences)):
        sequences[i] = sequences[i].replace('L', 'I')
        sequences[i] = sequences[i].replace('V', 'I')
        sequences[i] = sequences[i].replace('E', 'D')
        sequences[i] = sequences[i].replace('R', 'K')
    return sequences


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: {} <FASTA alignment>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1])
