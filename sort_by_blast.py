#!/usr/bin/env python3

import sys

import bioinfo

def read_blast(filename) :
    '''Read a tabular BLAST output file'''
    with open(filename) as infile :
        raw = infile.read()
        data = tuple(line.split('\t') for line in raw.split('\n') if len(line) > 0)
    return data

def match_blast_to_fasta(identifier, headers) :
    '''Get the index for the full header for the identifier from the BLAST file'''
    try :
        match = [header for header in headers if identifier in header][0]
    except IndexError :
        return False
    return headers.index(match)

def main() :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {} <BLAST tabular output> <FASTA file>\n'.format(sys.argv[0]))
        sys.exit(1)

    blast_data = read_blast(sys.argv[1])
    headers, seqs = bioinfo.read_fasta(sys.argv[2])
    
    for i in range(len(blast_data)) :
        ind = match_blast_to_fasta(blast_data[i][1], headers)
        print('>{}\n{}'.format(headers[ind], seqs[ind]))

if __name__ == '__main__' :
    main()
