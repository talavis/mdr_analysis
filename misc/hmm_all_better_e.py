#!/usr/bin/env python3

'''
Evaluates a HMMER output file based on an E value threshold.
Input: HMMER output file, E value threshold
Output: set of all accessions with E greater than the threshold
'''

import sys

def main(filename):
    with open(filename) as infile:
        current = ''
        identifier = ''
        first = True
        accs = set()
        for line in infile:
            if line[0] == '#':
                continue
            cols = line.split()
            if len(cols) < 2 :
                continue
            if cols[2] != current:
                current = cols[2]
                if identifier != current:
                    identifier = current
                    first = True
                    accs.add(identifier)
                if float(cols[4]) > float(sys.argv[2]) and first:
                    break
        print('\n'.join(list(accs)))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: <HMMER output file> <E threshold>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1])
