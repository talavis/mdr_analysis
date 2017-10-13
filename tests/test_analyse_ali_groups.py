#!/usr/bin/env python3
'''
Test analyse_ali_groups
'''

import analyse_ali_groups as aag


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
    assert aag.conservation(seqs) == expected


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile
    filename = tempfile.mkstemp()[1]
    data = ('>prot1\n' +
            'ACDEFGHIKLMNPQRSTVWY\n' +
            '>prot2\n' +
            'ACEEFAHIKIMNPSKSTVWY\n' +
            '>prot3\n' +
            'ACEEFGHIKVMNPQKSTVWA')
    with open(filename, 'w') as tmpfile:
        tmpfile.write(data)
    aag.main(filename)
    outerr = capsys.readouterr()
    expected = ('A\t1.00\nC\t1.00\nD\t1.00\nD\t1.00\nF\t1.00\n' +
                'G\t0.67\nH\t1.00\nI\t1.00\nK\t1.00\nI\t1.00\n' +
                'M\t1.00\nN\t1.00\nP\t1.00\nQ\t0.67\nK\t1.00\n' +
                'S\t1.00\nT\t1.00\nI\t1.00\nW\t1.00\nY\t0.67\n')

    assert outerr[0] == expected


def test_transform_prop():
    '''
    Test transform_prop()
    '''
    sequences = ['ACDEFGHIKLMNPQRSTVWY',
                 'ACDEFGHIKLMNPQRSTVWY',
                 'ACDEFGHIKLMNPQRSTVWY']
    expected = ['ACDDFGHIKIMNPQKSTIWY',
                'ACDDFGHIKIMNPQKSTIWY',
                'ACDDFGHIKIMNPQKSTIWY']

    assert aag.transform_prop(sequences) == expected
