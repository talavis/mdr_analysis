#!/usr/bin/env python3
'''Sort the sequences in a FASTA file according to
the order in a tabbed BLAST output file'''

import sys

import bioinfo


def read_blast(filename):
    '''Read a tabular BLAST output file'''
    with open(filename) as infile:
        raw = infile.read()
        data = tuple(line.split('\t') for line in raw.split('\n') if len(line) > 0)

    return data


def test_read_blast():
    '''Test read_blast'''
    import tempfile

    indata = '''sp|Q9BV79|MECR_HUMAN	gi|544346134|ref|NP_057095.3|	99.732	373	1	0	1	373	1	373	0.0	763
sp|Q9BV79|MECR_HUMAN	gi|397515849|ref|XP_003828155.1|	99.196	373	3	0	1	373	1	373	0.0	757'''
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    real = (['sp|Q9BV79|MECR_HUMAN', 'gi|544346134|ref|NP_057095.3|',
             '99.732', '373', '1', '0', '1', '373', '1', '373', '0.0', '763'],
            ['sp|Q9BV79|MECR_HUMAN', 'gi|397515849|ref|XP_003828155.1|',
             '99.196', '373', '3', '0', '1', '373', '1', '373', '0.0', '757'])

    assert read_blast(filename) == real


def match_blast_to_fasta(identifier, headers):
    '''Get the index for the full header for the identifier from the BLAST file.
       Returns False if not found'''
    match = [header for header in headers if identifier in header]
    if len(match) > 1:
        sys.stderr.write('I: identifier ({}) found in multiple headers; returning the first\n'.format(identifier))
    try:
        return headers.index(match[0])
    except IndexError:
        sys.stderr.write('I: identifier ({}) not found in any header\n'.format(identifier))
        return False


def test_match_blast_to_fasta():
    '''Test match_blast_to_fasta()'''
    headers = ('gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]',
               'gi|18408069|ref|NP_566881.1| putative trans-2-enoyl-CoA reductase [Arabidopsis thaliana]')
    assert match_blast_to_fasta('gi|18408069', headers) == 1
    assert match_blast_to_fasta('incorrect', headers) is False


def main(blast_file, filename):
    '''Sort the sequences in a FASTA file according to
    the order in a tabbed BLAST output file'''

    blast_data = read_blast(blast_file)
    headers, seqs = bioinfo.read_fasta(filename)

    for i in range(len(blast_data)):
        ind = match_blast_to_fasta(blast_data[i][1], headers)
        if ind is not False:
            print('>{}\n{}'.format(headers[ind], bioinfo.beautify_fasta(seqs[ind])))


def test_main(capsys):
    '''Test main()'''
    import tempfile

    indata_blast = '''sp|Q9BV79|MECR_HUMAN	gi|0000001|ref|NP_000001.1|	99.732	373	1	0	1	373	1	373	0.0	763
sp|Q9BV79|MECR_HUMAN	gi|0000003|ref|NP_000003.1|	99.196	373	3	0	1	373	1	373	0.0	757'''
    filename_blast = tempfile.mkstemp()[1]
    with open(filename_blast, 'w') as tmpf:
        tmpf.write(indata_blast)

    indata_fasta = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]
ACDEFGHIKL
>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]
ACDEAGHIKL
>gi|0000002|ref|NP_000002.1| Made-up data [Homo Sapiens]
ACDEFGHIKL
>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]
ADDEEGHILL'''

    filename_fasta = tempfile.mkstemp()[1]
    with open(filename_fasta, 'w') as tmpf:
        tmpf.write(indata_fasta)

    real = '''>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]
ACDEAGHIKL
>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]
ADDEEGHILL
'''
    main(filename_blast, filename_fasta)
    out, err = capsys.readouterr()
    assert out == real


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <BLAST tabular output> <FASTA file>\n'.format(sys.argv[0]))
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
