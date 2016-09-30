#!/usr/bin/env python3
'''
Tests for the bioinfo module
'''

import bioinfo


def test_beautify_fasta():
    '''Test beautify_fasta()'''
    seq = ('MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD' +
           'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD')
    assert bioinfo.beautify_fasta(seq) == ('MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD' + '\n' +
                                   'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD')
    assert bioinfo.beautify_fasta(seq[:15]) == 'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD'[:15]


def test_get_species():
    '''Test get_species(); two cases should work and one should fail'''
    ncbi_head = 'gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, mitochondrial precursor [Rattus norvegicus]'
    uniprot_head = 'sp|Q9BV79|MECR_HUMAN Trans-2-enoyl-CoA reductase, mitochondrial OS=Homo sapiens GN=MECR PE=1 SV=2'
    ensembl_head = 'MECR-001 peptide: ENSP00000263702 pep:KNOWN_protein_coding'
    assert bioinfo.get_species(ncbi_head) == 'Rattus norvegicus'
    assert bioinfo.get_species(uniprot_head) == 'Homo sapiens'
    assert bioinfo.get_species(ensembl_head) == 'Unknown'


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

    assert bioinfo.read_fasta(filename) == (realhead, realseq)