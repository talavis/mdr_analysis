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
    
def main(filename, REFSEQ) :
    headers, seqs = bioinfo.read_fasta(filename)

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
    filter_species(headers, seqs)

    for i in range(len(headers)) :
        print('>{}'.format(headers[i]))
        print(bioinfo.beautify_fasta(seqs[i]))

def test_main(capsys) :
    import tempfile

    indata = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]
ACDEFGHIKL
>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]
ACDEAGHIKL
>gi|0000002|ref|NP_000002.1| Made-up data [Homo Sapiens]
ACD
>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]
ACDEFGHIKL'''
    
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as f :
        f.write(indata)

    REAL = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]\nACDEFGHIKL\n>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]\nACDEAGHIKL\n'''
    main(filename, 'NP_000001.1')
    out, err = capsys.readouterr()
    assert out == REAL
        
if __name__ == '__main__' :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {} <seqfile> <refseq>\n'.format(sys.argv[0]))
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
