#!/usr/bin/env python3
'''Some commonly used bioinformatics funtions'''

import sys


def beautify_fasta(sequence, per_line=60):
    '''Limit the number of characters per line for a sequence'''
    i = per_line
    while i < len(sequence):
        sequence = sequence[:i] + '\n' + sequence[i:]
        i += per_line + 1

    return sequence


def test_beautify_fasta():
    '''Test beautify_fasta()'''
    seq = ('MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD' +
           'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD')
    assert beautify_fasta(seq) == ('MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD' + '\n' +
                                   'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD')
    assert beautify_fasta(seq[:15]) == 'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD'[:15]


def get_species(header):
    '''Get the species of a sequence'''
    # NCBI
    if header[:2] == 'gi':
        return header[header.index('[')+1:header.index(']')]
    # UniProt
    if header[:2] in ['sp', 'tr']:
        return header[header.index('OS=')+3:header.index('=', header.index('OS=')+3)-2].strip()

    sys.stderr.write('E: Unable to identify species({})'.format(header))

    return 'Unknown'


def test_get_species():
    '''Test get_species(); two cases should work and one should fail'''
    ncbi_head = 'gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]'
    uniprot_head = 'sp|Q9BV79|MECR_HUMAN Trans-2-enoyl-CoA reductase, mitochondrial OS=Homo sapiens GN=MECR PE=1 SV=2'
    ensembl_head = 'MECR-001 peptide: ENSP00000263702 pep:KNOWN_protein_coding'
    assert get_species(ncbi_head) == 'Rattus norvegicus'
    assert get_species(uniprot_head) == 'Homo sapiens'
    assert get_species(ensembl_head) == 'Unknown'


def read_fasta(filename):
    '''Read a FASTA file and return the headers and sequences as lists'''
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


def test_read_fasta():
    '''Test read_fasta() by reading a badly formatted FASTA file'''
    import tempfile

    indata = '''>gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]
MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQLKNLELTAVEG
    SDVHVKMLAAPINPSDINMIQGNYGLLPKLPAVGGNEGVGQVIAVGSSVSGLKPGDWVIPANAGLGTWRT
EAVFSEEALIGVPKDIPLQSAATLGVNPCTAYRMLVDFEQLQPGDSVIQNASNSGVGQAVIQIASALGLK
    TINVIRDRPDIKKLTDRLKDLGADYVLTEEELRMPETKNIFKDLPLPRLALNCVGGKSSTELLRHLAPGG
TMVTYGGMAKQPVTASVSMLIFKDLKLRGFWLSQWKKNHSPDEF
KELILILCNLIRQGQLTAPAWSGIPL
QDYQQALEASMKPFVSLKQILTM

>gi|18408069|ref|NP_566881.1| putative trans-2-enoyl-CoA reductase [Arabidopsis thaliana]
    MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVTRLVNLPPVEV
KENDVCVKMIAAPINPSDINRIEGVYPVRPPVPAVGGYEGVGEVYAVGSNVNGFSPGDWVIPSPPSSGTW


QTYVVKEESVWHKIDKECPMEYAATITVNPLTALRMLEDFVNLNSGDSVVQNGATSIVGQCVIQLARLRG
ISTINLIRDRAGSDEAREQLKALGADEVFSESQLNVKNVKSLLGNLPEPALGFNCVGGNAASLVLKYLRE
GGTMVTYGGMSKKPITVSTTSFIFKDLALRGFWLQSWLSMGKVKECREMIDYLLGLARDGKLKYETELVP
FEEFPVALDKALGKLGRQPKQVITF'''

    realhead = ['gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]',
                'gi|18408069|ref|NP_566881.1| putative trans-2-enoyl-CoA reductase [Arabidopsis thaliana]']
    realseq = ['MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQLKNLELTAVEGSDVHVKMLAAPINPSDINMIQGNYGLLPKLPAVGGNEGVGQVIAVGSSVSGLKPGDWVIPANAGLGTWRTEAVFSEEALIGVPKDIPLQSAATLGVNPCTAYRMLVDFEQLQPGDSVIQNASNSGVGQAVIQIASALGLKTINVIRDRPDIKKLTDRLKDLGADYVLTEEELRMPETKNIFKDLPLPRLALNCVGGKSSTELLRHLAPGGTMVTYGGMAKQPVTASVSMLIFKDLKLRGFWLSQWKKNHSPDEFKELILILCNLIRQGQLTAPAWSGIPLQDYQQALEASMKPFVSLKQILTM',
               'MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVTRLVNLPPVEVKENDVCVKMIAAPINPSDINRIEGVYPVRPPVPAVGGYEGVGEVYAVGSNVNGFSPGDWVIPSPPSSGTWQTYVVKEESVWHKIDKECPMEYAATITVNPLTALRMLEDFVNLNSGDSVVQNGATSIVGQCVIQLARLRGISTINLIRDRAGSDEAREQLKALGADEVFSESQLNVKNVKSLLGNLPEPALGFNCVGGNAASLVLKYLREGGTMVTYGGMSKKPITVSTTSFIFKDLALRGFWLQSWLSMGKVKECREMIDYLLGLARDGKLKYETELVPFEEFPVALDKALGKLGRQPKQVITF']

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    assert read_fasta(filename) == (realhead, realseq)


if __name__ == '__main__':
    sys.exit()
