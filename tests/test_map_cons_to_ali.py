#!/usr/bin/env python3
'''
Map conservation data to an alignment
'''

import map_cons_to_ali as mca


def test_map_cons(cons):
    '''
    Test map_cons()
    '''
    heads = ['prot1', 'prot2', 'prot3']
    seqs = ['ACDEF', 'ADDEF', 'ACDEH']
    ali_data = [0.0]*5
    data = ('prot3', (2, 3, 5), (0.3, 0.5, 0.7))
    expected = [0.0, 0.3, 0.5, 0.0, 0.7]
    assert mca.map_cons(heads, seqs, ali_data, data) == expected


def test_map_pos(cons):
    '''
    Test map_pos()
    '''
    assert False


def test_main():
    '''
    Test main()
    '''
    assert False


def test_read_data(capsys):
    '''
    Test read_data()
    '''
    import tempfile

    indata = ('# refprot\n' +
              '1\tA\tA\t0.7\n' +
              '2\tC\tC\t0.9\n' +
              '3\tD\tD\t0.95\n' +
              '4\tC\tE\t0.1\n' +
              '5\tE\tE\t1\n' +
              '6\t\tF\t0.5\n')


    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    positions = [1, 2, 3, 4, 5, 6]
    ratios = [0.7, 0.9, 0.95, 0.1, 1, 0.5]
    expected = ('refprot', positions, ratios)
    assert mca.read_data(filename) == expected
    
