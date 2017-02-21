#!/usr/bin/env python3
'''Calculate the conservation in an alignment'''

import sys

from Bio import AlignIO
from Bio.Align import AlignInfo


def get_most_conserved(freq_table, align_len):
    '''
    Determine the most conserved residue and its conservation
    rate in each position.
    Input should be a PSSM object or a list of dicts
    (PSSM.pssm, but without [0] in each tuple)
    '''
    result = list()
    # depending on input type
    try:
        num_positions = len(freq_table)
    except TypeError:
        num_positions = len(freq_table.pssm)

    for pos in range(num_positions):
        values = list(freq_table[pos].values())
        if values.count(max(values)) == 1:
            freq_table_inv = dict((j, i) for i, j in freq_table[pos].items())
            best_conserved = freq_table_inv[max(freq_table_inv)]
            score = freq_table[pos][best_conserved]/align_len
        else:
            best_conserved = 'X'
            score = max(values)/align_len
        result.append((score, best_conserved))

    return result


def make_freq_table(alignment):
    '''
    Make a pssm of an alignment
    '''
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus()
    freq_table = summary.pos_specific_score_matrix(consensus)

    return freq_table


def main(filename, refseq=None):
    '''
    Read an alignment in FASTA format
    Calculate the conservation per position
    '''
    # load alignments
    alignment = AlignIO.read(filename, 'fasta')

    headers = [s.name for s in alignment]
    if refseq is not None:
        try:
            refind = headers.index([h for h in headers if refseq in h][0])
        except IndexError:
            error = ('E: The reference sequence ({}) ' +
                     'not found among the sequences\n').format(refseq)
            sys.stderr.write(error)
            return False

    freq_table = make_freq_table(alignment)
    cons = get_most_conserved(freq_table, len(alignment))

    if refseq is None:
        refseq_p = ' '*len(alignment[0].seq)
    else:
        refseq_p = str(alignment[refind].seq)
        print('# {}'.format(refseq))
    pos = 1
    for i in range(len(freq_table.pssm)):
        if refseq_p[i] != '-':
            print('{ps}\t{rs}\t{mc}\t{rate:.3}'.format(ps=pos,
                                                       rs=refseq_p[i],
                                                       mc=cons[i][1],
                                                       rate=cons[i][0]))
            pos += 1


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.stderr.write(('Usage: {0} '.format(sys.argv[0]) +
                          '<alignment file> [reference seq]\n'))
        sys.exit(1)

    if len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        main(sys.argv[1])
