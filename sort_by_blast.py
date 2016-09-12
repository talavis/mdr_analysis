#!/usr/bin/env python3

import sys

import bioinfo

def read_blast(filename) :
    '''Read a tabular BLAST output file'''
    with open(filename) as infile :
        raw = infile.read()
        data = tuple(line.split('\t') for line in raw.split('\n') if len(line) > 0)
    return data
        
def main() :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {} <BLAST tabular output> <FASTA file>\n'.format(sys.argv[0]))
        sys.exit(1)

    blast_data = read_blast(sys.argv[1])
    seqs, headers = bioinfo.read_fasta(sys.argv[2])

if __name__ == '__main__' :
    main()
