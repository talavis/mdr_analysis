#!/usr/bin/env python3

import sys

import bioinfo

def read_blast(filename) :
    '''Read a tabular BLAST output file'''
    with open(filename) as infile :
        raw = infile.read()
        data = tuple(line.split('\t') for line in raw.split('\n') if len(line) > 0)
    return data

def test_read_blast() :
    import tempfile

    indata = '''sp|Q9BV79|MECR_HUMAN	gi|544346134|ref|NP_057095.3|	99.732	373	1	0	1	373	1	373	0.0	763
sp|Q9BV79|MECR_HUMAN	gi|397515849|ref|XP_003828155.1|	99.196	373	3	0	1	373	1	373	0.0	757
'''
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as f :
        f.write(indata)

    assert read_blast(filename) == (['sp|Q9BV79|MECR_HUMAN', 'gi|544346134|ref|NP_057095.3|',
                                     '99.732', '373', '1', '0', '1', '373', '1', '373', '0.0', '763'],
                                    ['sp|Q9BV79|MECR_HUMAN', 'gi|397515849|ref|XP_003828155.1|',
                                     '99.196', '373', '3', '0', '1', '373', '1', '373', '0.0', '757'])


def match_blast_to_fasta(identifier, headers) :
    '''Get the index for the full header for the identifier from the BLAST file.
       Returns False if not found'''
    try :
        match = [header for header in headers if identifier in header]
    except IndexError :
        sys.stderr.write('I: identifier ({}) not found in any header\n'.format(identifier))
        return False
    if len(match) > 1 :
        sys.stderr.write('I: identifier ({}) found in multiple headers; returning the first\n'.format(identifier))
    return headers.index(match[0])

def test_match_blast_to_fasta() :
    headers = ('gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]', 'gi|18408069|ref|NP_566881.1| putative trans-2-enoyl-CoA reductase [Arabidopsis thaliana]')
    assert match_blast_to_fasta('gi|18408069', headers) == 1

def main() :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {} <BLAST tabular output> <FASTA file>\n'.format(sys.argv[0]))
        sys.exit(1)

    blast_data = read_blast(sys.argv[1])
    headers, seqs = bioinfo.read_fasta(sys.argv[2])
    
    for i in range(len(blast_data)) :
        ind = match_blast_to_fasta(blast_data[i][1], headers)
        print('>{}\n{}'.format(headers[ind], bioinfo.beautify_fasta(seqs[ind])))

if __name__ == '__main__' :
    main()
