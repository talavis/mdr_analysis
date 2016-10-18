#!/usr/bin/env python3
'''Calculate the conservation in an alignment'''

import sys

import conservation

from Bio import AlignIO
from Bio.Align import AlignInfo


def helper_test_getalign():
    '''Generate a FASTA file for use in the testing functions'''
    import tempfile

    indata = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]
ACDEFGHIKL
>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]
ACDEAGHIKL
>gi|0000002|ref|NP_000002.1| Made-up data [Homo Sapiens]
ACDEFGHIKL
>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]
ADDEEGHILL'''

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    return filename


def test_get_most_conserved():
    '''Test get_most_conserved'''
    freq_table = [{'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 4.0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 3.0, 'F': 0, 'D': 1.0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 4.0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 4.0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 1.0,
                   'C': 0, 'F': 2.0, 'D': 0, 'G': 0, 'E': 1.0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 4.0, 'E': 0},
                  {'K': 0, 'H': 4.0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 4.0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 3.0, 'H': 0, 'I': 0, 'L': 1.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 4.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}]

    real = [(1.0, 'A'), (0.75, 'C'), (1.0, 'D'), (1.0, 'E'), (0.5, 'F'),
            (1.0, 'G'), (1.0, 'H'), (1.0, 'I'), (0.75, 'K'), (1.0, 'L')]

    assert conservation.get_most_conserved(freq_table, 4) == real


def test_make_freq_table():
    '''Test make_freq_table'''
    real = [('A', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 4.0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('C', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 3.0, 'F': 0, 'D': 1.0, 'G': 0, 'E': 0}),
            ('D', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 4.0, 'G': 0, 'E': 0}),
            ('E', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 4.0}),
            ('X', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 1.0,
                   'C': 0, 'F': 2.0, 'D': 0, 'G': 0, 'E': 1.0}),
            ('G', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 4.0, 'E': 0}),
            ('H', {'K': 0, 'H': 4.0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('I', {'K': 0, 'H': 0, 'I': 4.0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('K', {'K': 3.0, 'H': 0, 'I': 0, 'L': 1.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('L', {'K': 0, 'H': 0, 'I': 0, 'L': 4.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0})]

    alignment = AlignIO.read(helper_test_getalign(), 'fasta')

    assert conservation.make_freq_table(alignment).pssm == real


def test_main(capsys):
    '''Test main()'''
    real = ('A\tA\t1.0\nC\tC\t0.75\nD\tD\t1.0\nE\tE\t1.0\nA\tF\t0.5\n' +
            'G\tG\t1.0\nH\tH\t1.0\nI\tI\t1.0\nK\tK\t0.75\nL\tL\t1.0\n')
    conservation.main(helper_test_getalign(), 'NP_000001.1')
    out, err = capsys.readouterr()
    assert out == real

    real = (' \tA\t1.0\n \tC\t0.75\n \tD\t1.0\n \tE\t1.0\n \tF\t0.5\n' +
            ' \tG\t1.0\n \tH\t1.0\n \tI\t1.0\n \tK\t0.75\n \tL\t1.0\n')
    conservation.main(helper_test_getalign(), None)
    out, err = capsys.readouterr()
    assert out == real
    conservation.main(helper_test_getalign())
    out, err = capsys.readouterr()
    assert out == real
