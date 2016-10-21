#!/usr/bin/env python3
'''
Collect all positions with conservation greater than limit
Default: 0.9
'''

import sys


def main(filename, limit=0.9):
    '''
    Collect all positions with conservation greater than limit
    '''
    with open(filename) as infile:
        raw = infile.read()
    lines = raw.split('\n')
    if lines[0][0] == '#':
        refline = lines[0]
        lines = lines[1:]
    else:
        refline = None
    lines = [l for l in lines if len(l) > 0]
    hits = [l for l in lines
            if float(l.split('\t')[3]) > limit]

    if refline is not None:
        print(refline)
    for i in range(len(hits)):
        print('{}{}'.format(int(hits[i][0]), hits[i][1:]))


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.stderr.write('Usage: {} <conservation data> [limit=0.9]\n'.format(sys.argv[0]))
        sys.exit(1)

    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main(sys.argv[1], float(sys.argv[3]))
