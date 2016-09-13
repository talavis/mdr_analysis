#!/usr/bin/env python3

import sys

import bioinfo

def filter_length(headers, sequences, reflen, lmin = 0.9, lmax = 1.1) :
    '''In-place removal of all sequences that are outside the limits'''
    i = 0
    while i < len(sequences) :
        if reflen / len(sequences[i]) < lmin or reflen / len(sequences[i]) > lmax :
            headers.pop(i)
            sequences.pop(i)
        else :
            i += 1


def test_filter_length() :
    headers = ['header1', 'header2', 'header3', 'header4']
    sequences = ['ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW', 'FGHIKLMNPQRSTVWY', 'AAAAAGIRSTVWY']
    answer = (headers[:2], sequences[:2])
    filter_length(headers, sequences, len(sequences[0]))
    assert (headers, sequences) == answer

def filter_species(headers, sequences) :
    '''In-place removal of any extra sequences from a species. Assumes first hit is best hit.'''
    species = list()
    i = 0
    while i < len(sequences) :
        spec = bioinfo.get_species(headers[i])
        if spec in species :
            headers.pop(i)
            sequences.pop(i)
        else :
            species.append(spec)
            i += 1

def test_filter_species() :
    headers = ['gi|header1 [Homo sapiens]', 'gi|header2 [Homo sapiens]',
               'sp|header3 OS=Rattus norvegicus GN=Unknown', 'gi|header4 [Rattus norvegicus]']
    sequences = ['ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVW', 'FGHIKLMNPQRSTVWY', 'AAAAAGIKLMNPQRSTVWY']
    answer = ([headers[0], headers[2]], [sequences[0], sequences[2]])
    filter_species(headers, sequences)
    assert (headers, sequences) == answer
    
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
    filter_length(headers, seqs, REFLEN)    
    
if __name__ == '__main__' :
    main()
