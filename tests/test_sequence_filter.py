#!/usr/bin/env python3
'''
Tests for the sequence_filter module
'''

import sys

import sequence_filter
import bioinfo


def test_filter_nonsense():
    '''Test filter_nonsense'''
    headers = ['header1', 'header2', 'header3', 'header4']
    sequences = ['ACDEXXXXXXXXPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW',
                 'FGXXXXQRSTVWYACDEHKLMN', 'AAAAAGIRSTVWY']
    answer = (headers[1:], sequences[1:])
    sequence_filter.filter_nonsense(headers, sequences)
    assert (headers, sequences) == answer
    headers = ['header1', 'header2', 'header3', 'header4']
    sequences = ['ACDEXXXXXXXXPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW',
                 'FGXXXXQRSTVWYACDEHKLMN', 'AAAAAGIRSTVWY']
    answer = ([headers[1], headers[3]], [sequences[1], sequences[3]])
    sequence_filter.filter_nonsense(headers, sequences, 0.0)
    assert (headers, sequences) == answer
    headers = ['header1', 'header2', 'header3', 'header4']
    sequences = ['ACDEXXXXXXXXPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW',
                 'FGXXXXQRSTVWYACDEHKLMN', 'AAAAAGIRSTVWY']
    answer = (headers[:], sequences[:])
    sequence_filter.filter_nonsense(headers, sequences, 1.0)
    assert (headers, sequences) == answer


def test_filter_length():
    '''Test filter_length'''
    headers = ['header1', 'header2', 'header3', 'header4']
    sequences = ['ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW', 'FGHIKLMNPQRSTVWY', 'AAAAAY']
    answer = (headers[:3], sequences[:3])
    sequence_filter.filter_length(headers, sequences, len(sequences[0]))
    assert (headers, sequences) == answer


def test_filter_species():
    '''Test filter_species()'''
    headers = ['gi|header1 [Homo sapiens]', 'gi|header2 [Homo sapiens]',
               'sp|header3 OS=Rattus norvegicus GN=Unknown', 'gi|header4 [Rattus norvegicus]']
    sequences = ['ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW',
                 'FGHIKLMNPQRSTVWY', 'AAAAAGIKLMNPQRSTVWY']
    answer = ([headers[0], headers[2]], [sequences[0], sequences[2]])
    sequence_filter.filter_species(headers, sequences)
    assert (headers, sequences) == answer


def test_main(capsys):
    '''Testing main()'''
    import tempfile

    indata = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]
ACDEFGHIKL
>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]
ACXXXXXIKL
>gi|0000002|ref|NP_000002.1| Made-up data [Homo Sapiens]
ACD
>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]
ACDEFGHIKL'''

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    real = ('>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]\n'
            + 'ACDEFGHIKL\n'
            + '>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]\n'
            + 'ACDEFGHIKL\n')

    sequence_filter.main(filename, 'NP_000001.1')
    out, err = capsys.readouterr()
    assert out == real
