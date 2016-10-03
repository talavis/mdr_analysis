#!/usr/bin/env python3

'''
Identify the first non-matching hit
Input: HMMER output file
Output: The first protein giving a different best hit
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
                if cols[0].lower() != identifier and first:
                    print('{}\t{}'.format(current, cols[4]))
                    first = False
    
        
if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: {} <HMMER output file>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1])
