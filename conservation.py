#!/usr/bin/env python3
'''Trim an alignment to only show columns with aa from specific class of ADH'''

import sys

from Bio import AlignIO
from Bio.Align import AlignInfo

def make_freq_table(alignment) :
    '''Make a pssm of an alignment'''
    summary1 = AlignInfo.SummaryInfo(alignment)
    consensus = summary1.dumb_consensus()
    freq_table = summary1.pos_specific_score_matrix(consensus)

    return freq_table

def calculate_conservation(freq_table, align_len) :
    '''Determine the most conserved residue and its conservation rate in each position'''
    result = list()
    for p in range(len(freq_table.pssm)) :
        freq_table_inv = dict((j, i) for i, j in freq_table[p].items())
        best_conserved = freq_table_inv[max(freq_table_inv)]
        score = freq_table[p][best_conserved]/align_len
        result.append((score, best_conserved))

    return result

def test_make_pssm() :
    pass

def main() :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {0} <Alignment> <reference sequence>\n'.format(sys.argv[0]))
        sys.exit(1)

    REFSEQ = sys.argv[2]

    # load alignments
    filename = sys.argv[1]
    alignment = AlignIO.read(filename, 'fasta')

    headers = [s.name for s in alignment]
    REFIND = headers.index([h for h in headers if REFSEQ in h][0])
    SEQLEN = len(str(alignment[REFIND].seq).replace('-', ''))

    freq_table = make_freq_table(alignment)
    cons = calculate_conservation(freq_table, len(alignment))

    for i in range(len(freq_table.pssm)) :
        print(freq_table[i])
        print(cons[i])

if __name__ == '__main__' :
    main()
