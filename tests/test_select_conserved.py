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

    indata = ('A\tA\t0.7\n' +
              'C\tC\t0.9\n' +
              'D\tD\t0.95\n' +
              'C\tE\t0.1\n' +
              'E\tE\t1\n' +
              'F\tF\t0.5\n')

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
