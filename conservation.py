#!/usr/bin/env python3
'''Calculate the conservation in an alignment'''

import sys

from Bio import AlignIO
from Bio.Align import AlignInfo


def get_most_conserved(freq_table, align_len):
    '''Determine the most conserved residue and its conservation rate in each position
    Input should be a PSSM object or a list of dicts (PSSM.pssm, but without [0] in each tuple)'''
    result = list()
    # depending on input type
    try:
        num_positions = len(freq_table)
    except TypeError:
        num_positions = len(freq_table.pssm)

    for pos in range(num_positions):
        freq_table_inv = dict((j, i) for i, j in freq_table[pos].items())
        best_conserved = freq_table_inv[max(freq_table_inv)]
        score = freq_table[pos][best_conserved]/align_len
        result.append((score, best_conserved))

    return result


def make_freq_table(alignment):
    '''Make a pssm of an alignment'''
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus()
    freq_table = summary.pos_specific_score_matrix(consensus)

    return freq_table


def main(filename, refseq):
    '''Read an alignment in FASTA format
    Calculate the conservation per position'''
    # load alignments
    alignment = AlignIO.read(filename, 'fasta')

    headers = [s.name for s in alignment]
    try:
        refind = headers.index([h for h in headers if refseq in h][0])
    except IndexError:
        sys.stderr.write('E: The reference sequence ({}) not found among the sequences\n'.format(refseq))
        sys.exit(1)

    seqlen = len(str(alignment[refind].seq).replace('-', ''))

    freq_table = make_freq_table(alignment)
    cons = get_most_conserved(freq_table, len(alignment))

    for i in range(len(freq_table.pssm)):
        print(cons[i])


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {0} <alignment file> <reference sequence>\n'.format(sys.argv[0]))
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
