#!/usr/bin/env python3
'''
Map conservation data to an alignment
'''

import map_cons_to_ali as mca


def test_map_cons():
    '''
    Test map_cons()
    '''
    heads = ['prot1', 'prot2', 'prot3']
    seqs = ['ACDEF', 'ADDEF', 'ACDEH']
    data = ('prot3', (2, 3, 5), (0.3, 0.5, 0.7))
    expected = [0.0, 0.3, 0.5, 0.0, 0.7]
    assert mca.map_cons(heads, seqs, data) == expected
    seqs = ['AC-DEF', 'AD-DEF', 'ACDE-H']
    data = ('prot2', (2, 3, 5), (0.3, 0.5, 0.7))
    expected = [0.0, 0.3, 0.0, 0.5, 0.0, 0.7]
    assert mca.map_cons(heads, seqs, data) == expected


def test_map_pos(capsys):
    '''
    Test map_pos()
    '''
    pos = 7
    seq = '----ACFER----FR'
    assert mca.map_pos(pos, seq) == 14
    pos = 1
    seq = '----ACFER----FR'
    assert mca.map_pos(pos, seq) == 4
    pos = 10
    assert mca.map_pos(pos, seq) is False
    out, err = capsys.readouterr()
    assert err == 'E: position 10 not found in ----ACFER----FR\n'


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    # correct
    indata1 = ('# prot2\n' +
               '1\tA\tA\t0.7\n' +
               '2\tC\tC\t0.9\n' +
               '3\tD\tD\t0.95\n' +
               '4\tC\tE\t0.1\n' +
               '5\tE\tE\t1\n' +
               '6\t\tF\t0.5\n')
    data1_name = tempfile.mkstemp()[1]
    with open(data1_name, 'w') as tmpf:
        tmpf.write(indata1)

    indata2 = ('# prot1\n' +
               '1\tA\tA\t0.1\n' +
               '2\tC\tC\t0.5\n' +
               '3\tD\tD\t0.82\n' +
               '4\tC\tE\t0.13\n' +
               '5\tE\tE\t1.0\n' +
               '6\t\tF\t0.05\n')
    data2_name = tempfile.mkstemp()[1]
    with open(data2_name, 'w') as tmpf:
        tmpf.write(indata2)

    inseq = ('>prot1\n' +
             'ACCDEF-\n' +
             '>prot2\n' +
             'A-DDEFG\n' +
             '>prot3\n' +
             'A--DEFG\n')
    fasta_name = tempfile.mkstemp()[1]
    with open(fasta_name, 'w') as tmpf:
        tmpf.write(inseq)

    mca.main(fasta_name, [data1_name, data2_name])
    out, err = capsys.readouterr()
    assert out == ('ali:prot1\tA\tC\tC\tD\tE\tF\t-\n' +
                   'ali:prot2\tA\t-\tD\tD\tE\tF\tG\n' +
                   'ali:prot3\tA\t-\t-\tD\tE\tF\tG\n' +
                   'cons:prot2\t0.7\t0.0\t0.9\t0.95\t0.1\t1.0\t0.5\n' +
                   'cons:prot1\t0.1\t0.5\t0.82\t0.13\t1.0\t0.05\t0.0\n')
    # no reference in data file
    indata = ('1\tA\tA\t0.1\n' +
              '2\tC\tC\t0.5\n' +
              '3\tD\tD\t0.82\n' +
              '4\tC\tE\t0.13\n' +
              '5\tE\tE\t1.0\n' +
              '6\t\tF\t0.05\n')
    data_name = tempfile.mkstemp()[1]
    with open(data_name, 'w') as tmpf:
        tmpf.write(indata)

    assert mca.main(fasta_name, [data_name]) is False


def test_read_data():
    '''
    Test read_data()
    '''
    import tempfile

    # correct
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

    # no reference in data file
    indata = ('1\tA\tA\t0.1\n' +
              '2\tC\tC\t0.5\n' +
              '3\tD\tD\t0.82\n' +
              '4\tC\tE\t0.13\n' +
              '5\tE\tE\t1.0\n' +
              '6\t\tF\t0.05\n')
    data_name = tempfile.mkstemp()[1]
    with open(data_name, 'w') as tmpf:
        tmpf.write(indata)

    assert mca.read_data(data_name) is False
