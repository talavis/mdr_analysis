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


def test_conservation():
    '''
    Test analyse_ali()
    '''
    seqs = ('ACDEF',
            'ACEEF',
            'ACDEE',
            'CCDDD',
            'DCADF')
    expected = list(zip('ACDEF', (0.6, 1.0, 0.6, 0.6, 0.6)))
    assert conservation(seqs) == expected



def main(filename):
    '''
    Read a FASTA alignment and perform a conservation analysis
    using groups of residues
    '''


def test_main(capsys):
    '''
    Test main()
    '''


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: {} <FASTA alignment>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1])
