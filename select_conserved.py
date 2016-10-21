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
    hits = [l for l in enumerate(lines)
            if len(l[1]) > 0 and
            float(l[1].split('\t')[2]) > limit]

    if refline is not None:
        print(refline)
    for i in range(len(hits)):
        print('{}\t{}'.format(hits[i][0] + 1, hits[i][1]))


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.stderr.write('Usage: {} <conservation data> [limit=0.9]\n'.format(sys.argv[0]))
        sys.exit(1)

    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main(sys.argv[1], float(sys.argv[3]))
