#!/usr/bin/env python3
'''
Filter the sequences, keeping only the ones present in 
the provided taxonomy file.
'''

import sys

import bioinfo


def main(fasta_file, taxonomy_file):
    '''
    Filter the sequences, keeping only the ones present in 
    the provided taxonomy file.
    '''
    species = read_species(taxonomy_file)
    headers, sequences = bioinfo.read_fasta(fasta_file)

    for i in range(len(headers)):
        if bioinfo.get_species(headers[i]) in species:
            print('>{}\n{}\n'.format(headers[i],
                                     sequences[i]))


def read_species(taxonomy_file):
    '''
    Read a UniProt taxonomy file
    Return a list of species names
    '''
    raw = open(taxonomy_file).read().split('\n')[1:]

    return [line.split('\t')[2]
            for line
            in raw
            if len(line) > 0]


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <fasta file> <taxonomy file>\n'.format(sys.argv[0]))
        sys.exit(1)
    if main(sys.argv[1], sys.argv[2]) is False:
        sys.exit(1)
