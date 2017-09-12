#!/usr/bin/env python3
'''
Testing conservation.py
'''

from Bio import AlignIO

import conservation

def helper_test_getalign():
    '''
    Generate a FASTA file for use in the testing functions
    '''
    import tempfile

    indata = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]
ACD--GHIKL
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
    '''
    Test get_most_conserved()
    '''
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

    freq_table = [{'-': 0, 'A': 4.0, 'C': 0, 'D': 0, 'E': 0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 0, 'A': 0, 'C': 3.0, 'D': 1.0, 'E': 0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 4.0, 'A': 0, 'C': 0, 'D': 0.0, 'E': 0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 1.0, 'A': 0, 'C': 0, 'D': 0, 'E': 3.0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 1.0, 'A': 1.0, 'C': 0, 'D': 0, 'E': 1.0,
                   'F': 1.0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 4.0, 'H': 0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 4.0, 'I': 0, 'K': 0, 'L': 0},
                  {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 0, 'I': 4.0, 'K': 0, 'L': 0},
                  {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 0, 'I': 0, 'K': 3.0, 'L': 1.0},
                  {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 4.0}]

    real = [(1.0, 'A'), (0.75, 'C'), (1.0, '-'), (0.75, 'E'), (0.25, 'X'),
            (1.0, 'G'), (1.0, 'H'), (1.0, 'I'), (0.75, 'K'), (1.0, 'L')]

    assert conservation.get_most_conserved(freq_table, 4) == real


def test_group_res_prop():
    '''
    Test group_res_prop()
    '''
    sequences = ['ACDEFGHIKLMNPQRSTVWY',
                 'ACDEFGHIKLMNPQRSTVWY',
                 'ACDEFGHIKLMNPQRSTVWY']
    expected = ['ACDDFGHIKIMNPNKSSIWY',
                'ACDDFGHIKIMNPNKSSIWY',
                'ACDDFGHIKIMNPNKSSIWY']

    assert conservation.group_res_prop(sequences) == expected


def test_main(capsys):
    '''
    Test main()
    '''
    # with a reference sequence
    real = ('# NP_000000.1\n' +
            '1\tA\tA\t1.0\n' +
            '2\tC\tC\t0.75\n' +
            '3\tD\tD\t1.0\n' +
            '4\tG\tG\t1.0\n' +
            '5\tH\tH\t1.0\n' +
            '6\tI\tI\t1.0\n' +
            '7\tK\tK\t0.75\n' +
            '8\tL\tL\t1.0\n')
    conservation.main(helper_test_getalign(), 'NP_000000.1')
    out, err = capsys.readouterr()
    assert out == real

    # Without reference sequence
    real = ('1\t \tA\t1.0\n' +
            '2\t \tC\t0.75\n' +
            '3\t \tD\t1.0\n' +
            '4\t \tE\t0.75\n' +
            '5\t \tX\t0.25\n' +
            '6\t \tG\t1.0\n' +
            '7\t \tH\t1.0\n' +
            '8\t \tI\t1.0\n' +
            '9\t \tK\t0.75\n' +
            '10\t \tL\t1.0\n')
    conservation.main(helper_test_getalign(), None)
    out, err = capsys.readouterr()
    assert out == real
    conservation.main(helper_test_getalign())
    out, err = capsys.readouterr()
    assert out == real

    # incorrect reference
    assert conservation.main(helper_test_getalign(), 'incorrect') is False
    out, err = capsys.readouterr()
    expected = 'E: The reference sequence (incorrect) not found among the sequences\n'
    assert err == expected

    import tempfile

    indata = ('>gi|0000000|ref|NP_000000.1| Made-up data ' +
              '[Rattus norvegicus]\n' +
              'ACD--GHIKL\n'
              '>gi|0000001|ref|NP_000001.1| Made-up data ' +
              '[Arabidopsis thaliana]\n' +
              'ACDEAGHIKL\n' +
              '>gi|0000002|ref|NP_000002.1| Made-up data ' +
              '[Homo Sapiens]\n' +
              'ACD--GHIKL\n' +
              '>gi|0000003|ref|NP_000003.1| Made-up data ' +
              '[Arabidopsis thaliana]\n' +
              'ADD--GHILL\n')

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)
    expected = ('# NP_000001.1\n' +
                '1\tA\tA\t1.0\n' +
                '2\tC\tC\t0.75\n' +
                '3\tD\tD\t1.0\n' +
                '4\tE\t-\t0.75\n' +
                '5\tA\t-\t0.75\n' +
                '6\tG\tG\t1.0\n' +
                '7\tH\tH\t1.0\n' +
                '8\tI\tI\t1.0\n' +
                '9\tK\tK\t0.75\n' +
                '10\tL\tL\t1.0\n')
    conservation.main(filename, 'NP_000001.1')
    out, err = capsys.readouterr()
    assert out == expected


def test_make_freq_table():
    '''
    Test make_freq_table()
    '''
    real = [('A', {'-': 0, 'A': 4.0, 'C': 0, 'D': 0, 'E': 0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0}),
            ('C', {'-': 0, 'A': 0, 'C': 3.0, 'D': 1.0, 'E': 0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0}),
            ('D', {'-': 0, 'A': 0, 'C': 0, 'D': 4.0, 'E': 0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0}),
            ('E', {'-': 1.0, 'A': 0, 'C': 0, 'D': 0, 'E': 3.0,
                   'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0}),
            ('X', {'-': 1.0, 'A': 1.0, 'C': 0, 'D': 0, 'E': 1.0,
                   'F': 1.0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0}),
            ('G', {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 4.0, 'H': 0, 'I': 0, 'K': 0, 'L': 0}),
            ('H', {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 4.0, 'I': 0, 'K': 0, 'L': 0}),
            ('I', {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 0, 'I': 4.0, 'K': 0, 'L': 0}),
            ('K', {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 0, 'I': 0, 'K': 3.0, 'L': 1.0}),
            ('L', {'-': 0, 'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0,
                   'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 4.0})]

    alignment = AlignIO.read(helper_test_getalign(), 'fasta')

    assert conservation.make_freq_table(alignment).pssm == real
