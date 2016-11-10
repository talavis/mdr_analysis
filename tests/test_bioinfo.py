#!/usr/bin/env python3
'''
Tests for the bioinfo module
'''

import bioinfo


def test_beautify_fasta():
    '''
    Test beautify_fasta()
    '''
    seq = ('MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD' +
           'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD')

    expected = ('MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD\n' +
                'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD')
    assert bioinfo.beautify_fasta(seq) == expected

    expected = 'MCTTGQVIRCKAAILWKPGAPFSIEEVEVAPPKAKEVRIKVVATGLCGTEMKVLGSKHLD'[:15]
    assert bioinfo.beautify_fasta(seq[:15]) == expected


def test_get_species():
    '''
    Test get_species(); two cases should work and one should fail
    '''
    ncbi_head = ('gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, ' +
                 'mitochondrial precursor [Rattus norvegicus]')
    uniprot_head = ('sp|Q9BV79|MECR_HUMAN Trans-2-enoyl-CoA reductase, ' +
                    'mitochondrial OS=Homo sapiens GN=MECR PE=1 SV=2')
    ensembl_head = 'MECR-001 peptide: ENSP00000263702 pep:KNOWN_protein_coding'
    assert bioinfo.get_species(ncbi_head) == 'Rattus norvegicus'
    assert bioinfo.get_species(uniprot_head) == 'Homo sapiens'
    assert bioinfo.get_species(ensembl_head) == 'Unknown'


def test_get_structseq(capsys):
    '''
    Test get_structseq
    '''
    e_heads = ['1U3W:A|PDBID|CHAIN|SEQUENCE',
               '1U3W:B|PDBID|CHAIN|SEQUENCE']
    e_seqs = [('STAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKM' +
               'VAAGICRSDEHVVSGNLVTPLPVILGHEAAGIVESVGEGV' +
               'TTVKPGDKVIPLFTPQCGKCRICKNPESNYCLKNDLGNPR' +
               'GTLQDGTRRFTCSGKPIHHFVGVSTFSQYTVVDENAVAKI' +
               'DAASPLEKVCLIGCGFSTGYGSAVKVAKVTPGSTCAVFGL' +
               'GGVGLSVVMGCKAAGAARIIAVDINKDKFAKAKELGATEC' +
               'INPQDYKKPIQEVLKEMTDGGVDFSFEVIGQLDTMMASLL' +
               'CCHEACGTSVIVGVPPDSQNLSINPMLLLTGRTWKGAIFG' +
               'GFKSKESVPKLVADFMAKKFSLDALITNVLPFEKINEGFD' +
               'LLRSGKSIRTVLTF'),
              ('STAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKM' +
               'VAAGICRSDEHVVSGNLVTPLPVILGHEAAGIVESVGEGV' +
               'TTVKPGDKVIPLFTPQCGKCRICKNPESNYCLKNDLGNPR' +
               'GTLQDGTRRFTCSGKPIHHFVGVSTFSQYTVVDENAVAKI' +
               'DAASPLEKVCLIGCGFSTGYGSAVKVAKVTPGSTCAVFGL' +
               'GGVGLSVVMGCKAAGAARIIAVDINKDKFAKAKELGATEC' +
               'INPQDYKKPIQEVLKEMTDGGVDFSFEVIGQLDTMMASLL' +
               'CCHEACGTSVIVGVPPDSQNLSINPMLLLTGRTWKGAIFG' +
               'GFKSKESVPKLVADFMAKKFSLDALITNVLPFEKINEGFD' +
               'LLRSGKSIRTVLTF')]

    assert bioinfo.get_structseq('1u3w') == (e_heads, e_seqs)
    assert bioinfo.get_structseq('abcde') is False
    err = capsys.readouterr()[1]
    assert err == 'E: could not retrieve sequence for structure abcde\n'


def test_read_fasta():
    '''
    Test read_fasta() by reading a badly formatted FASTA file
    '''
    import tempfile

    indata = ('>gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, ' +
              'mitochondrial precursor [Rattus norvegicus]\n' +
              'MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQLKNLELTAVEG\n' +
              '    SDVHVKMLAAPINPSDINMIQGNYGLLPKLPAVGGNEGVGQVIAVGSSVSGLKPGDWVIPANAGLGTWRT\n' +
              'EAVFSEEALIGVPKDIPLQSAATLGVNPCTAYRMLVDFEQLQPGDSVIQNASNSGVGQAVIQIASALGLK\n' +
              '    TINVIRDRPDIKKLTDRLKDLGADYVLTEEELRMPETKNIFKDLPLPRLALNCVGGKSSTELLRHLAPGG\n' +
              'TMVTYGGMAKQPVTASVSMLIFKDLKLRGFWLSQWKKNHSPDEF\n' +
              'KELILILCNLIRQGQLTAPAWSGIPL\n' +
              'QDYQQALEASMKPFVSLKQILTM\n' +
              ' \n' +
              '>gi|18408069|ref|NP_566881.1| putative ' +
              'trans-2-enoyl-CoA reductase [Arabidopsis thaliana]\n' +
              '    MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVTRLVNLPPVEV\n' +
              'KENDVCVKMIAAPINPSDINRIEGVYPVRPPVPAVGGYEGVGEVYAVGSNVNGFSPGDWVIPSPPSSGTW\n' +
              ' \n' +
              '\n' +
              'QTYVVKEESVWHKIDKECPMEYAATITVNPLTALRMLEDFVNLNSGDSVVQNGATSIVGQCVIQLARLRG\n' +
              'ISTINLIRDRAGSDEAREQLKALGADEVFSESQLNVKNVKSLLGNLPEPALGFNCVGGNAASLVLKYLRE\n' +
              'GGTMVTYGGMSKKPITVSTTSFIFKDLALRGFWLQSWLSMGKVKECREMIDYLLGLARDGKLKYETELVP\n' +
              'FEEFPVALDKALGKLGRQPKQVITF\n')

    realhead = [('gi|8393848|ref|NP_058905.1| trans-2-enoyl-CoA reductase, ' +
                 'mitochondrial precursor [Rattus norvegicus]'),
                ('gi|18408069|ref|NP_566881.1| putative trans-2-enoyl-CoA ' +
                 'reductase [Arabidopsis thaliana]')]
    realseq = [('MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQL' +
                'KNLELTAVEGSDVHVKMLAAPINPSDINMIQGNYGLLPKLPAVGGNEGVGQVIAVGSSVS' +
                'GLKPGDWVIPANAGLGTWRTEAVFSEEALIGVPKDIPLQSAATLGVNPCTAYRMLVDFEQ' +
                'LQPGDSVIQNASNSGVGQAVIQIASALGLKTINVIRDRPDIKKLTDRLKDLGADYVLTEE' +
                'ELRMPETKNIFKDLPLPRLALNCVGGKSSTELLRHLAPGGTMVTYGGMAKQPVTASVSML' +
                'IFKDLKLRGFWLSQWKKNHSPDEFKELILILCNLIRQGQLTAPAWSGIPLQDYQQALEAS' +
                'MKPFVSLKQILTM'),
               ('MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVT' +
                'RLVNLPPVEVKENDVCVKMIAAPINPSDINRIEGVYPVRPPVPAVGGYEGVGEVYAVGSN' +
                'VNGFSPGDWVIPSPPSSGTWQTYVVKEESVWHKIDKECPMEYAATITVNPLTALRMLEDF' +
                'VNLNSGDSVVQNGATSIVGQCVIQLARLRGISTINLIRDRAGSDEAREQLKALGADEVFS' +
                'ESQLNVKNVKSLLGNLPEPALGFNCVGGNAASLVLKYLREGGTMVTYGGMSKKPITVSTT' +
                'SFIFKDLALRGFWLQSWLSMGKVKECREMIDYLLGLARDGKLKYETELVPFEEFPVALDK' +
                'ALGKLGRQPKQVITF')]

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    assert bioinfo.read_fasta(filename) == (realhead, realseq)


def test_read_fasta_raw():
    '''
    Tests read_fasta_raw()
    '''
    indata = ('>1U3W:A|PDBID|CHAIN|SEQUENCE\n' +
              'STAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKM' +
              'VAAGICRSDEHVVSGNLVTPLPVILGHEAAGIVESVGEGV\n' +
              'TTVKPGDKVIPLFTPQCGKCRICKNPESNYCLKNDLGNPR' +
              'GTLQDGTRRFTCSGKPIHHFVGVSTFSQYTVVDENAVAKI\n' +
              'DAASPLEKVCLIGCGFSTGYGSAVKVAKVTPGSTCAVFGL' +
              'GGVGLSVVMGCKAAGAARIIAVDINKDKFAKAKELGATEC\n' +
              'INPQDYKKPIQEVLKEMTDGGVDFSFEVIGQLDTMMASLL' +
              'CCHEACGTSVIVGVPPDSQNLSINPMLLLTGRTWKGAIFG\n' +
              'GFKSKESVPKLVADFMAKKFSLDALITNVLPFEKINEGFD' +
              'LLRSGKSIRTVLTF\n' +
              '>1U3W:B|PDBID|CHAIN|SEQUENCE\n' +
              'STAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKM' +
              'VAAGICRSDEHVVSGNLVTPLPVILGHEAAGIVESVGEGV\n' +
              'TTVKPGDKVIPLFTPQCGKCRICKNPESNYCLKNDLGNPR' +
              'GTLQDGTRRFTCSGKPIHHFVGVSTFSQYTVVDENAVAKI\n' +
              'DAASPLEKVCLIGCGFSTGYGSAVKVAKVTPGSTCAVFGL' +
              'GGVGLSVVMGCKAAGAARIIAVDINKDKFAKAKELGATEC\n' +
              'INPQDYKKPIQEVLKEMTDGGVDFSFEVIGQLDTMMASLL' +
              'CCHEACGTSVIVGVPPDSQNLSINPMLLLTGRTWKGAIFG\n' +
              'GFKSKESVPKLVADFMAKKFSLDALITNVLPFEKINEGFD' +
              'LLRSGKSIRTVLTF\n')
    e_heads = ['1U3W:A|PDBID|CHAIN|SEQUENCE',
               '1U3W:B|PDBID|CHAIN|SEQUENCE']
    e_seqs = [('STAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKM' +
               'VAAGICRSDEHVVSGNLVTPLPVILGHEAAGIVESVGEGV' +
               'TTVKPGDKVIPLFTPQCGKCRICKNPESNYCLKNDLGNPR' +
               'GTLQDGTRRFTCSGKPIHHFVGVSTFSQYTVVDENAVAKI' +
               'DAASPLEKVCLIGCGFSTGYGSAVKVAKVTPGSTCAVFGL' +
               'GGVGLSVVMGCKAAGAARIIAVDINKDKFAKAKELGATEC' +
               'INPQDYKKPIQEVLKEMTDGGVDFSFEVIGQLDTMMASLL' +
               'CCHEACGTSVIVGVPPDSQNLSINPMLLLTGRTWKGAIFG' +
               'GFKSKESVPKLVADFMAKKFSLDALITNVLPFEKINEGFD' +
               'LLRSGKSIRTVLTF'),
              ('STAGKVIKCKAAVLWELKKPFSIEEVEVAPPKAHEVRIKM' +
               'VAAGICRSDEHVVSGNLVTPLPVILGHEAAGIVESVGEGV' +
               'TTVKPGDKVIPLFTPQCGKCRICKNPESNYCLKNDLGNPR' +
               'GTLQDGTRRFTCSGKPIHHFVGVSTFSQYTVVDENAVAKI' +
               'DAASPLEKVCLIGCGFSTGYGSAVKVAKVTPGSTCAVFGL' +
               'GGVGLSVVMGCKAAGAARIIAVDINKDKFAKAKELGATEC' +
               'INPQDYKKPIQEVLKEMTDGGVDFSFEVIGQLDTMMASLL' +
               'CCHEACGTSVIVGVPPDSQNLSINPMLLLTGRTWKGAIFG' +
               'GFKSKESVPKLVADFMAKKFSLDALITNVLPFEKINEGFD' +
               'LLRSGKSIRTVLTF')]

    assert bioinfo.read_fasta_raw(indata) == (e_heads, e_seqs)


def test_res_at_pos(capsys):
    '''
    Test res_at_pos()
    '''
    # correct
    seq = 'ACDEFGHIKLMNPQRSTVWY'
    pos = (0, 4, 8, 10, 15, 19)
    expected = ['A', 'F', 'K', 'M', 'S', 'Y']
    assert bioinfo.res_at_pos(seq, pos) == expected
    # incorrect position, too high
    pos = (0, 4, 8, 10, 15, 25)
    assert bioinfo.res_at_pos(seq, pos) is False
    err = capsys.readouterr()[1]
    e_err = 'E: requested positions outside the protein sequence\n'
    assert err == e_err
    # incorrect position, too low
    pos = (-1, 4, 8, 10, 15, 19)
    assert bioinfo.res_at_pos(seq, pos) is False
    err = capsys.readouterr()[1]
    assert err == e_err
