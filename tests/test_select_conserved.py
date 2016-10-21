#!/usr/bin/env python3
'''
Collect all positions with conservation greater than limit
Default: 0.9
'''

import select_conserved as sc


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    # without reference
    indata = ('1\tA\tA\t0.7\n' +
              '2\tC\tC\t0.9\n' +
              '3\tD\tD\t0.95\n' +
              '4\tC\tE\t0.1\n' +
              '5\tE\tE\t1\n' +
              '6\tF\tF\t0.5\n')

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    # data read
    sc.main(filename)
    out, err = capsys.readouterr()
    assert out == ('3\tD\tD\t0.95\n' +
                   '5\tE\tE\t1\n')
    # limit
    sc.main(filename, 0.6)
    out, err = capsys.readouterr()
    assert out == ('1\tA\tA\t0.7\n' +
                   '2\tC\tC\t0.9\n' +
                   '3\tD\tD\t0.95\n' +
                   '5\tE\tE\t1\n')

    # with reference sequence
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

    # data read
    sc.main(filename)
    out, err = capsys.readouterr()
    assert out == ('# refprot\n' +
                   '3\tD\tD\t0.95\n' +
                   '5\tE\tE\t1\n')
