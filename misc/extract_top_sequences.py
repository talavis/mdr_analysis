#!/usr/bin/env python3

'''
Extract the top X sequences from a FASTA file
'''

import sys

import bioinfo


def main(filename, nr_sequences):
    headers, sequences = bioinfo.read_fasta(filename)

    for i in range(int(nr_sequences)):
        print('>{}'.format(headers[i]))
        print(bioinfo.beautify_fasta(sequences[i]))

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <FASTA file> <#sequences>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

        
