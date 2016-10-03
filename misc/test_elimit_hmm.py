#!/usr/bin/env python3

'''
Evaluates a HMMER output file based on an E value threshold.
Input: HMMER output file, E value threshold
Output: List of the last matching sequence for each of the 
HMM models and the last hits' E value
'''

import sys

def main(filename):
    with open(filename) as infile:
        current = ''
        identifier = ''
        first = True
        for line in infile:
            if line[0] == '#':
                continue
            cols = line.split()
            if len(cols) < 2 :
                continue
            if cols[2] != current:
                current = cols[2]
                if identifier != current[:current.index('_')]:
                    identifier = current[:current.index('_')]
                    first = True
                if float(cols[4]) > float(sys.argv[2]) and first:
                    print('{}\t{}'.format(current, cols[4]))
                    first = False

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: <HMMER output file> <E threshold>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1])
