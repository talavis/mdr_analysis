#!/usr/bin/env python3
'''
Filter the BLAST results based on the last hit with the correct gene name

Assume that all sequence headers are formatted UniProt style (GN=)
Asssume that the first sequence has the correct gene name

Returns all sequences above the last sequence with correct gene name
'''

import sys

import bioinfo


def main(hits_fasta_file):
    '''
    Filter the BLAST results based on the last hit with the correct gene name
    '''
    headers, seqs = bioinfo.read_fasta(hits_fasta_file)

    wanted = headers[0][headers[0].index('GN=')+3:
                        headers[0].index(' ', headers[0].index('GN='))]
    hits = [h for h in headers if wanted in h]

    limit = headers.index(hits[-1])
    for i in range(limit+1):
        print('>' + headers[i])
        print(bioinfo.beautify_fasta(seqs[i]))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: {} <BLAST hits.fasta>\n'.format(sys.argv[0]))
        sys.exit(1)
    if main(sys.argv[1]) is False:
        sys.exit(1)
