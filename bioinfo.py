#!/usr/bin/env python3
'''
Some commonly used bioinformatics funtions
'''

import sys


def beautify_fasta(sequence, per_line=60):
    '''
    Limit the number of characters per line for a sequence
    '''
    i = per_line
    while i < len(sequence):
        sequence = sequence[:i] + '\n' + sequence[i:]
        i += per_line + 1

    return sequence


def get_species(header):
    '''
    Get the species of a sequence
    '''
    # NCBI
    if header[:2] == 'gi':
        return header[header.index('[')+1:header.index(']')]
    # UniProt
    if header[:2] in ['sp', 'tr']:
        return header[header.index('OS=')+3:
                      header.index('=', header.index('OS=')+3)-2].strip()

    sys.stderr.write('E: Unable to identify species({})'.format(header))

    return 'Unknown'


def read_fasta(filename):
    '''
    Read a FASTA file.
    Return the headers and sequences as lists
    '''
    with open(filename) as infile:
        seqs = list()
        headers = list()
        for line in infile:
            if line[0] == '>':
                headers.append(line[1:].strip())
                seqs.append('')
            elif len(line.strip()) > 0:
                seqs[-1] += line.strip()

    return (headers, seqs)


def read_fasta_raw(raw):
    '''
    Parse a FASTA formatted text string.
    Return the headers and sequences as lists
    '''
    seqs = list()
    headers = list()
    for line in raw.split('\n'):
        if len(line) == 0:
            continue
        if line[0] == '>':
            headers.append(line[1:].strip())
            seqs.append('')
        elif len(line.strip()) > 0:
            seqs[-1] += line.strip()

    return (headers, seqs)
