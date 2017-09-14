#!/usr/bin/env python3
'''
Convert a contact table from ICM to an ICM selection
'''

import sys


def cont_to_sele(filename):
    '''
    Convert a contact table from ICM to an ICM selection
    '''
    with open(filename) as infile:
        data = infile.read().split('\n')
    seles = [row.split(',')[0] for row in data if row and row[0] != 'A']
    base = seles[0][:seles[0].index('.')]
    sele = 'a_{}.a/{}'.format(base, ','.join([sele[sele.index('^'):]
                                            for sele in seles]))
    print(sele)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: {} <ICM contacts file>\n'.format(sys.argv[0]))
        sys.exit(1)
    cont_to_sele(sys.argv[1])

