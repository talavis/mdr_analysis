#!/usr/bin/env python3
'''
Retrieve the wanted accession codes from a postree
'''

import re
import sys


def find_accs(data, pattern, revert=False):
    '''
    Find accs in a "postree"
    '''
    data = data.split('|')
    # uniprot accession code
    regex = re.compile('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}')
    # should be in the format X_ACC
    if revert:
        matches = [dat for dat in data if pattern not in dat]
    else:
        matches = [dat for dat in data if pattern in dat]
    accs = [acc[2:] for acc in matches if regex.fullmatch(acc[2:])]
    return accs


def test_find_accs():
    '''
    Test find_accs()
    '''
    data = ('(tr|L_Q0G839|L_Q0G839_CAEEL:0.1,tr|M_A0A091MYM1|M_A0A091MYM1_9PASS:0.2,' +
            '(tr|L_H2M6E8|L_H2M6E8_ORYLA:0.3,tr|L_C3Z2I6|L_C3Z2I6_BRAFL:0.4)' +
            '1.00:0.5)0.897;')
    expected = ['Q0G839', 'H2M6E8', 'C3Z2I6']
    assert find_accs(data, 'L_') == expected
    expected = ['A0A091MYM1']
    assert find_accs(data, 'L_', True) == expected


def main(filename, pattern, revert=False):
    '''
    Read a file with "postree" and find the accession codes
    with the residue specifed in pattern
    '''
    data = open(filename).read()
    accs = find_accs(data, pattern, revert)

    for acc in accs:
        print(acc)


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile
    filename = tempfile.mkstemp()[1]
    data = ('(tr|L_Q0G839|L_Q0G839_CAEEL:0.1,tr|L_A0A091MYM1|L_A0A091MYM1_9PASS:0.2,' +
            '(tr|L_H2M6E8|L_H2M6E8_ORYLA:0.3,tr|M_C3Z2I6|M_C3Z2I6_BRAFL:0.4)' +
            '1.00:0.5)0.897;')
    with open(filename, 'w') as tmpfile:
        tmpfile.write(data)
    main(filename, 'L_')
    expected = 'Q0G839\nA0A091MYM1\nH2M6E8\n'
    assert capsys.readouterr()[0] == expected
    main(filename, 'M_')
    expected = 'C3Z2I6\n'
    assert capsys.readouterr()[0] == expected
    main(filename, 'C_')
    expected = ''
    assert capsys.readouterr()[0] == expected
    # revert
    main(filename, 'L_', True)
    expected = 'C3Z2I6\n'


if __name__ == '__main__':
    if len(sys.argv) not in (3, 4):
        sys.stderr.write('Usage: {} <tree file> <pattern> [revert]\n'.format(sys.argv[0]))
        sys.exit(1)

    if len(sys.argv) == 4 and sys.argv[3] == 'revert':
        main(sys.argv[1], sys.argv[2], True)
    else:
        main(sys.argv[1], sys.argv[2])
