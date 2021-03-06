#!/usr/bin/env python3
'''Filter a list of sequences according to different criteria'''

import sys

import bioinfo


def filter_nonsense(headers, sequences, unk_rate=0.1):
    '''
    In-place removal of all sequences that contain a lot of nonsense, e.g. X
    '''
    i = 0
    while i < len(sequences):
        if sequences[i].count('X') / len(sequences[i].replace('-', '')) > unk_rate:
            headers.pop(i)
            sequences.pop(i)
        else:
            i += 1


def filter_length(headers, sequences, reflen, lmin=0.5, lmax=1.5):
    '''In-place removal of all sequences that are outside the supplied limits'''
    i = 0
    while i < len(sequences):
        if reflen / len(sequences[i]) < lmin or reflen / len(sequences[i]) > lmax:
            headers.pop(i)
            sequences.pop(i)
        else:
            i += 1


def filter_species(headers, sequences):
    '''In-place removal of any extra sequences from a species. Assumes first hit is best hit.'''
    species = list()
    i = 0
    while i < len(sequences):
        spec = bioinfo.get_species(headers[i])
        if spec in species:
            headers.pop(i)
            sequences.pop(i)
        else:
            species.append(spec)
            i += 1


def main(filename, refseq):
    '''Read a FASTA file with sequences and remove sequences
    that do not pass the filters'''
    headers, seqs = bioinfo.read_fasta(filename)

    # confirm that there is _one_ matching sequence
    refseq_matches = [s for s in headers if refseq in s]
    if not refseq_matches:
        error = 'E: Reference sequence ({}) not found in the sequence file\n'.format(refseq)
        sys.stderr.write(error)
        return False
    if len(refseq_matches) > 1:
        error = ('E: Reference sequence ({}) matches '.format(refseq) +
                 'multiple sequences in the sequence file\n' +
                 'E: Match example: {}\n'.format(refseq_matches[0]) +
                 'E: Match example: {}\n'.format(refseq_matches[1]))
        sys.stderr.write(error)
        return False

    reflen = len(seqs[headers.index(refseq_matches[0])])
    filter_length(headers, seqs, reflen)
    filter_nonsense(headers, seqs)
    filter_species(headers, seqs)


    for i in range(len(headers)):
        print('>{}'.format(headers[i]))
        print(bioinfo.beautify_fasta(seqs[i]))


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <seqfile> <refseq>\n'.format(sys.argv[0]))
        sys.exit(1)

    if main(sys.argv[1], sys.argv[2]) is False:
        sys.exit(1)
