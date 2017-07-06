#!/usr/bin/env python3
'''Sort the sequences in a FASTA file according to
the order in a tabbed BLAST output file'''

import sys

import bioinfo


def read_blast(filename):
    '''Read a tabular BLAST output file'''
    with open(filename) as infile:
        raw = infile.read()
        data = tuple(line.split('\t') for line in raw.split('\n') if line)

    return data


def match_blast_to_fasta(identifier, headers):
    '''Get the index for the full header for the identifier from the BLAST file.
       Returns False if not found'''
    match = [header for header in headers if identifier in header]
    if len(match) > 1:
        sys.stderr.write('I: identifier ({}) found in '.format(identifier) +
                         'multiple headers; returning the first\n')
    try:
        return headers.index(match[0])
    except IndexError:
        sys.stderr.write('I: identifier ({}) not found in any header\n'.format(identifier))
        return False


def main(blast_file, filename):
    '''Sort the sequences in a FASTA file according to
    the order in a tabbed BLAST output file'''

    blast_data = read_blast(blast_file)
    headers, seqs = bioinfo.read_fasta(filename)

    for i in range(len(blast_data)):
        ind = match_blast_to_fasta(blast_data[i][1], headers)
        if ind is not False:
            print('>{}\n{}'.format(headers[ind], bioinfo.beautify_fasta(seqs[ind])))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <BLAST tabular output> <FASTA file>\n'.format(sys.argv[0]))
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
