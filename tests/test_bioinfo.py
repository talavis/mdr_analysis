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


def test_parse_icm_sel():
    '''
    Test parse_icm_sel()
    '''
    indata = 'a_1yb5.b/^Q167,^R170,^L190:^G193,^W294,^K296,^V298\n'

    expected = [167, 170, 190, 191, 192, 193, 294, 296, 298]

    assert bioinfo.parse_icm_sel(indata) == expected

    # incorrect code
    indata = 'a_1yb5asd.b/^Q167,^R170,^L190:^G193,^W294,^K296,^V298\n'
    assert bioinfo.parse_icm_sel(indata) is False


def test_read_fasta(capsys):
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

    # incorrect file name
    assert bioinfo.read_fasta('__incorrect') is False
    err = capsys.readouterr()[1]
    assert err == 'E: Cannot find FASTA file: __incorrect\n'


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


def test_read_icm_res():
    '''
    Test read_icm_res()
    '''
    import tempfile

    indata = ('a_1yb5.3,13/*|a_1yb5.a/^R257,^M260:^A261|' +
              'a_1yb5.b/^R12,^H36,^V50:^E51,^I54,^Y59:^R61,' +
              '^P63,^Y67:^S71,^S85:^F87,^I131,^S274,^K276,^Q287\n' +
              'a_1yb5.b/^Q167,^R170,^L190:^G193,^W294,^K296,^V298\n' +
              'a_1yb5.a/^K210,^V213:^E215,^G217,^S235:^L237,^K262\n')

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    expected = [[12, 36, 50, 51, 54, 59, 60, 61, 63, 67, 68, 69, 70,
                 71, 85, 86, 87, 131, 257, 260, 261, 274, 276, 287],
                [167, 170, 190, 191, 192, 193, 294, 296, 298],
                [210, 213, 214, 215, 217, 235, 236, 237, 262]]

    assert bioinfo.read_icm_res(filename) == expected


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
