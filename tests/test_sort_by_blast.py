#!/usr/bin/env python3
'''
Tests for the sort_by_blast module
'''

import bioinfo
import sort_by_blast


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

    assert sort_by_blast.read_blast(filename) == real


def test_match_blast_to_fasta():
    '''Test match_blast_to_fasta()'''
    headers = ('gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]',
               'gi|18408069|ref|NP_566881.1| putative trans-2-enoyl-CoA reductase [Arabidopsis thaliana]')
    assert sort_by_blast.match_blast_to_fasta('gi|18408069', headers) == 1
    assert sort_by_blast.match_blast_to_fasta('incorrect', headers) is False


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
    sort_by_blast.main(filename_blast, filename_fasta)
    out, err = capsys.readouterr()
    assert out == real
