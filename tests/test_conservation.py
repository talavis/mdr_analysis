#!/usr/bin/env python3
'''
Testing conservation.py
'''

from Bio import AlignIO

import conservation
import pytest


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


def test_conservation(capsys):
    '''
    Test conservation()
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
    conservation.conservation(helper_test_getalign(), 'NP_000000.1')
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
    conservation.conservation(helper_test_getalign(), refseq=None)
    out, err = capsys.readouterr()
    assert out == real
    conservation.conservation(helper_test_getalign())
    out, err = capsys.readouterr()
    assert out == real

    # incorrect reference
    assert conservation.conservation(helper_test_getalign(), 'incorrect') is False
    err = capsys.readouterr()[1]
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
    conservation.conservation(filename, 'NP_000001.1')
    out = capsys.readouterr()[0]
    assert out == expected

    # with residue grouping
    real = ('1\t \tA\t1.0\n' +
            '2\t \tC\t0.75\n' +
            '3\t \tD\t1.0\n' +
            '4\t \tD\t0.75\n' +
            '5\t \tX\t0.25\n' +
            '6\t \tG\t1.0\n' +
            '7\t \tH\t1.0\n' +
            '8\t \tI\t1.0\n' +
            '9\t \tK\t0.75\n' +
            '10\t \tI\t1.0\n')
    conservation.conservation(helper_test_getalign(), group_res=True)
    out = capsys.readouterr()[0]
    assert out == real

    # with ignore gaps
    real = ('1\t \tA\t1.0\n' +
            '2\t \tC\t0.75\n' +
            '3\t \tD\t1.0\n' +
            '4\t \tE\t1.0\n' +
            '5\t \tX\t0.333\n' +
            '6\t \tG\t1.0\n' +
            '7\t \tH\t1.0\n' +
            '8\t \tI\t1.0\n' +
            '9\t \tK\t0.75\n' +
            '10\t \tL\t1.0\n')
    conservation.conservation(helper_test_getalign(), ign_gaps=True)
    out = capsys.readouterr()[0]
    assert out == real


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

    assert conservation.get_most_conserved(freq_table, ign_gaps=True) == real

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

    assert conservation.get_most_conserved(freq_table) == real
    # with ign_gaps=True
    real = [(1.0, 'A'), (0.75, 'C'), (0.0, 'X'), (1.0, 'E'), (1/3, 'X'),
            (1.0, 'G'), (1.0, 'H'), (1.0, 'I'), (0.75, 'K'), (1.0, 'L')]

    assert conservation.get_most_conserved(freq_table, ign_gaps=True) == real


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


def test_parse_parameters():
    '''
    Test parse_parameters()
    '''
    # correct parameters
    try:
        params = ('filename.ali')
        expected = ('filename.ali', None, False, False)
        assert conservation.parse_parameters(params) == expected
    except ValueError:
        print('error')
    params = ('filename.ali', 'reference=P123456')
    expected = ('filename.ali', 'P123456', False, False)
    assert conservation.parse_parameters(params) == expected
    params = ('filename.ali', 'reference=P123456', 'group_res=y')
    expected = ('filename.ali', 'P123456', True, False)
    assert conservation.parse_parameters(params) == expected
    params = ('filename.ali', 'group_res=y')
    expected = ('filename.ali', None, True, False)
    assert conservation.parse_parameters(params) == expected
    params = ('filename.ali', 'group_res=y', 'ign_gaps=Y')
    expected = ('filename.ali', None, True, True)
    assert conservation.parse_parameters(params) == expected
    # different order
    params = ('filename.ali', 'group_res=y', 'ign_gaps=y', 'reference=P123456')
    expected = ('filename.ali', 'P123456', True, True)
    assert conservation.parse_parameters(params) == expected
    # incorrect parameters
    with pytest.raises(ValueError) as excinfo:
        expected = 'E: incorrect parameter (group_res=a)'
        params = ('filename.ali', 'group_res=a')
        conservation.parse_parameters(params)
    assert str(excinfo.value) == expected
    with pytest.raises(ValueError) as excinfo:
        expected = 'E: incorrect parameter (ign_gaps=s)'
        params = ('filename.ali', 'ign_gaps=s')
        conservation.parse_parameters(params)
    assert str(excinfo.value) == expected
    with pytest.raises(ValueError) as excinfo:
        expected = 'E: incorrect parameter (abc)'
        params = ('filename.ali', 'abc')
        conservation.parse_parameters(params)
    assert str(excinfo.value) == expected


def test_print_use(capsys):
    '''
    Test print_use()
    '''
    expected = ('Usage: conservation.py <alignment file> ' +
                '[reference=seq] [group_res=y/N] [ign_gaps=y/N]\n')
    conservation.print_use('conservation.py')
    err = capsys.readouterr()[1]
    assert err == expected
