#!/usr/bin/env python3

import sys

import bioinfo

def filter_length(headers, sequences, reflen, lmin = 0.9, lmax = 1.1) :
    '''Remove all sequences that are outside the limits'''
    i = 0
    while i < len(sequences) :
        if reflen / len(sequences[i]) < lmin or reflen / len(sequences[i]) > lmax :
            headers.pop(i)
            sequences.pop(i)
        i += 1
    return headers, sequences

def main() :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {} <seqfile> <refseq>\n'.format(sys.argv[0]))
        sys.exit(1)

    REFSEQ = sys.argv[2]
    headers, seqs = bioinfo.read_fasta(sys.argv[1])

    # confirm that there is _one_ matching sequence
    refseq_matches = [s for s in headers if REFSEQ in s]
    if len(refseq_matches) == 0 :
        sys.stderr.write('E: Reference sequence ({}) not found in the sequence file\n'.format(REFSEQ))
        sys.exit(1)
    if len(refseq_matches) > 1 :
        sys.stderr.write('E: Reference sequence ({}) matches multiple sequences in the sequence file\n'.format(REFSEQ))
        sys.stderr.write('E: Match example: {}\n'.format(refseq_matches[0]))
        sys.stderr.write('E: Match example: {}\n'.format(refseq_matches[1]))
        sys.exit(1)

    REFLEN = len(seqs[headers.index(refseq_matches[0])])
    headers, seqs = filter_length(headers, seqs, REFLEN)

    
    
if __name__ == '__main__' :
    main()
