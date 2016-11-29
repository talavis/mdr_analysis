#!/usr/bin/env python3
'''
Map active site data to an alignment
'''

import map_as_to_ali as maa


def test_find_refprot_index(capsys):
    '''
    Test find_refprot_index
    '''
    heads = ['1U3W_A|PDBID|CHAIN|SEQUENCE',
             '1U3W_B|PDBID|CHAIN|SEQUENCE',
             'PDBID|1U3W_C|CHAIN|SEQUENCE',
             '1U3W_D|PDBID|CHAIN|SEQUENCE']
    refprot = 'W_C'

    # success
    assert maa.find_refprot_index(refprot, heads) == 2
    # failure
    assert maa.find_refprot_index('asdf', heads) is False
    err = capsys.readouterr()[1]
    assert err == 'E: reference sequence asdf not found\n'


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    inasdata = ('- Num  Res. Type ---- SS Molecule ---- Object - sf - sfRatio\n' +
                '   1  ser   Amino    S _  a            pdbc     3.3  0.02  .  a_pdbc.a/^S1\n' +
                '   7  ile   Amino    I H  a            pdbc    37.0  0.16  .  a_pdbc.a/^I7\n')
    asdata_name = tempfile.mkstemp()[1]
    with open(asdata_name, 'w') as tmpf:
        tmpf.write(inasdata)

    # with fasta
    inseq = ('>prot1\n' +
             'STA-GKVIKCKAAVLW\n' +
             '>prot2\n' +
             'STAAGK-IKCKAAV-W\n' +
             '>prot3\n' +
             'S-ATGKLIKCKAAVL-\n')
    fasta_name = tempfile.mkstemp()[1]
    with open(fasta_name, 'w') as tmpf:
        tmpf.write(inseq)

    maa.main(fasta_name, 'prot1', '1u3w', asdata_name)
    out = capsys.readouterr()[0]
    assert out == ('ali:prot1\tS\tT\tA\t-\tG\tK\tV\tI\tK\tC\tK\tA\tA\tV\tL\tW\n' +
                   'ali:prot2\tS\tT\tA\tA\tG\tK\t-\tI\tK\tC\tK\tA\tA\tV\t-\tW\n' +
                   'ali:prot3\tS\t-\tA\tT\tG\tK\tL\tI\tK\tC\tK\tA\tA\tV\tL\t-\n' +
                   'as:prot1_1u3w\tX\t\t\t\t\t\t\tX\t\t\t\t\t\t\t\t\n')

    # with map
    inmap = ('ali:prot1\tS\tT\tA\t-\tG\tK\tV\tI\tK\tC\tK\tA\tA\tV\tL\tW\n' +
             'ali:prot2\tS\tT\tA\tA\tG\tK\t-\tI\tK\tC\tK\tA\tA\tV\t-\tW\n' +
             'ali:prot3\tS\t-\tA\tT\tG\tK\tL\tI\tK\tC\tK\tA\tA\tV\tL\t-\n' +
             'cons:prot2\t0.7\t0.0\t0.9\t0.95\t0.1\t1.0\t0.5\n' +
             'cons:prot1\t0.1\t0.5\t0.82\t0.13\t1.0\t0.05\t0.0\n')

    map_name = tempfile.mkstemp()[1]
    with open(map_name, 'w') as tmpf:
        tmpf.write(inmap)

    maa.main(map_name, 'prot1', '1u3w', asdata_name)
    out = capsys.readouterr()[0]
    assert out == ('ali:prot1\tS\tT\tA\t-\tG\tK\tV\tI\tK\tC\tK\tA\tA\tV\tL\tW\n' +
                   'ali:prot2\tS\tT\tA\tA\tG\tK\t-\tI\tK\tC\tK\tA\tA\tV\t-\tW\n' +
                   'ali:prot3\tS\t-\tA\tT\tG\tK\tL\tI\tK\tC\tK\tA\tA\tV\tL\t-\n' +
                   'cons:prot2\t0.7\t0.0\t0.9\t0.95\t0.1\t1.0\t0.5\n' +
                   'cons:prot1\t0.1\t0.5\t0.82\t0.13\t1.0\t0.05\t0.0\n' +
                   'as:prot1_1u3w\tX\t\t\t\t\t\t\tX\t\t\t\t\t\t\t\t\n')

    # incorrect structure
    assert maa.main(map_name, 'prot1', 'asdfg', asdata_name) is False
    capsys.readouterr()

    # incorrect reference
    assert maa.main(map_name, '1torp', '1u3w', asdata_name) is False
    capsys.readouterr()

    # incorrect icmfile
    assert maa.main(map_name, 'prot1', '1u3w', '______.txt') is False
    capsys.readouterr()

    # incorrect positions
    inasdata = ('- Num  Res. Type ---- SS Molecule ---- Object - sf - sfRatio\n' +
                '   1  ser   Amino    S _  a            pdbc     3.3  0.02  .  a_pdbc.a/^S1\n' +
                ' 100  ile   Amino    I H  a            pdbc    37.0  0.16  .  a_pdbc.a/^I100\n')
    asdata_name = tempfile.mkstemp()[1]
    with open(asdata_name, 'w') as tmpf:
        tmpf.write(inasdata)

    assert maa.main(map_name, 'prot1', '1u3w', asdata_name) is False


def test_map_as():
    '''
    Test map_as()
    '''
    # without gaps
    heads = ['prot1', 'prot2', 'prot3']
    seqs = ['ACDEF', 'ADDEF', 'ACDEH']
    ref = 'prot3'
    as_pos = [0, 1, 4]
    expected = ['X', 'X', '', '', 'X']
    assert maa.map_as(heads, seqs, ref, as_pos) == expected
    # with gaps
    heads = ['prot1', 'prot2', 'prot3']
    seqs = ['ACDEF-HIKLMNP-RST',
            '-CDEFGHIKLMNPQ-ST',
            'ACDE--HIKLMNPQRST']
    ref = 'prot3'
    as_pos = [1, 5, 10]
    expected = ['', 'X', '', '', '',
                '', '', 'X', '', '',
                '', '', 'X', '', '',
                '', '']
    assert maa.map_as(heads, seqs, ref, as_pos) == expected
    # incorrect reference
    assert maa.map_as(heads, seqs, 'incorrect', as_pos) is False


def test_map_pos(capsys):
    '''
    Test map_pos()
    '''
    pos = 7
    seq = '----ACFER----FR'
    assert maa.map_pos(pos, seq) == 14
    pos = 1
    seq = '----ACFER----FR'
    assert maa.map_pos(pos, seq) == 4
    pos = 10
    assert maa.map_pos(pos, seq) is False
    err = capsys.readouterr()[1]
    assert err == 'E: position 10 not found in ----ACFER----FR\n'


def test_map_sequences(capsys):
    '''
    Test map_sequences()
    '''
    struct_seq = ('MHHHHHHSSGVDLGTENLYFQSMMQKLVVTRLSPNFREAV' +
                  'TLSRDCPVPLPGDGDLLVRNRFVGVNASDINYSAGRYDPS' +
                  'VKPPFDIGFEGIGEVVALGLSASARYTVGQAVAYMAPGSF' +
                  'AEYTVVPASIATPVPSVKPEYLTLLVSGTTAYISLKELGG' +
                  'LSEGKKVLVTAAAGGTGQFAMQLSKKAKCHVIGTCSSDEK' +
                  'SAFLKSLGCDRPINYKTEPVGTVLKQEYPEGVDVVYESVG' +
                  'GAMFDLAVDALATKGRLIVIGFISGYQTPTGLSPVKAGTL' +
                  'PAKLLKKSASVQGFFLNHYLSKYQAAMSHLLEMCVSGDLV' +
                  'CEVDLGDLSPEGRFTGLESIFRAVNYMYMGKNTGKIVVEL' +
                  'PH')
    prot_seq = ('MLRLVPTGARAIVDMSYARHFLDFQGSAIP' +
                'QAMQKLVVTRLSPNFREAVTLSRDCPVPLP' +
                'GDGDLLVRNRFVGVNASDINYSAGRYDPSV' +
                'KPPFDIGFEGIGEVVALGLSASARYTVGQA' +
                'VAYMAPGSFAEYTVVPASIATPVPSVKPEY' +
                'LTLLVSGTTAYISLKELGGLSEGKKVLVTA' +
                'AAGGTGQFAMQLSKKAKCHVIGTCSSDEKS' +
                'AFLKSLGCDRPINYKTEPVGTVLKQEYPEG' +
                'VDVVYESVGGAMFDLAVDALATKGRLIVIG' +
                'FISGYQTPTGLSPVKAGTLPAKLLKKSASV' +
                'QGFFLNHYLSKYQAAMSHLLEMCVSGDLVC' +
                'EVDLGDLSPEGRFTGLESIFRAVNYMYMGK' +
                'NTGKIVVELPHSVNSKL')
    positions = [65, 71, 88, 105, 270, 284]
    residues = ['N', 'Y', 'F', 'Y', 'G', 'L']

    # correct mapping
    expected = [74, 80, 97, 114, 279, 293]
    assert maa.map_sequences(struct_seq, prot_seq, positions, residues) == expected
    capsys.readouterr()

    # incorrect residue
    positions = [65, 71, 105, 270, 284]
    residues = ['N', 'Y', 'A', 'G', 'L']
    assert maa.map_sequences(struct_seq, prot_seq, positions, residues) is False
    err = capsys.readouterr()[1]
    assert err == ('E: the protein sequence does not match the position data; ' +
                   'Position 106 should be A, but is Y\n')

    # incorrect mapping
    residues = ['H', 'Y', 'Y', 'G', 'L']
    positions = [2, 65, 71, 105, 270, 284]
    assert maa.map_sequences(struct_seq, prot_seq, positions, residues) is False
    err = capsys.readouterr()[1]
    assert err == 'E: the sequence HHHHHS is not found in the second protein\n'

    # position outside protein
    residues = ['N', 'Y', 'Y', 'G', 'L']
    positions = [65, 71, 105, 270, 500]
    assert maa.map_sequences(struct_seq, prot_seq, positions, residues) is False
    err = capsys.readouterr()[1]
    assert err == 'E: position outside protein (500)\n'


def test_read_icmdata(capsys):
    '''
    Test read_icmdata()
    '''
    import tempfile

    # correct file
    indata = ('- Num  Res. Type ---- SS Molecule ---- Object - sf - sfRatio\n' +
              '   66  asn   Amino    N _  a            pdbc     3.3  0.02  .  a_pdbc.a/^N66\n' +
              '   72  tyr   Amino    Y H  a            pdbc    37.0  0.16  .  a_pdbc.a/^Y72\n' +
              '  106  tyr   Amino    Y _  a            pdbc    20.5  0.09  .  a_pdbc.a/^Y106\n' +
              '  271  gly   Amino    G _  a            pdbc     0.8  0.01  .  a_pdbc.a/^G271\n' +
              '  285  leu   Amino    L H  b            pdbc    18.5  0.09  .  a_pdbc.b/^L285\n')
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    positions = [66, 72, 106, 271, 285]
    residues = ['N', 'Y', 'Y', 'G', 'L']
    expected = [positions, residues]
    assert maa.read_icmdata(filename) == expected

    # without header
    indata = ('   66  asn   Amino    N _  a            pdbc     3.3  0.02  .  a_pdbc.a/^N66\n' +
              '   72  tyr   Amino    Y H  a            pdbc    37.0  0.16  .  a_pdbc.a/^Y72\n' +
              '  106  tyr   Amino    Y _  a            pdbc    20.5  0.09  .  a_pdbc.a/^Y106\n' +
              '  271  gly   Amino    G _  a            pdbc     0.8  0.01  .  a_pdbc.a/^G271\n' +
              '  285  leu   Amino    L H  b            pdbc    18.5  0.09  .  a_pdbc.b/^L285\n')
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)
    assert maa.read_icmdata(filename) == expected

    # file not existing
    assert maa.read_icmdata('______.txt') is False
    err = capsys.readouterr()[1]
    assert err == 'E: file ______.txt not found\n'

    # incorrect formatting in file
    indata = ('- Num  Res. Type ---- SS Molecule ---- Object - sf - sfRatio\n' +
              '   6  6  asn   Amino    N _  a            pdbc     3.3  0.02  .  a_pdbc.a/^N66\n' +
              '   7 2  tyr   Amino    Y H  a            pdbc    37.0  0.16  .  a_pdbc.a/^Y72\n' +
              '  106  tyr   Amino    Y _  a            pdbc    20.5  0.09  .  a_pdbc.a/^Y106\n' +
              '271  gly   Amino    G _  a            pdbc     0.8  0.01  .  a_pdbc.a/^G271\n' +
              'incorrect')
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    assert maa.read_icmdata(filename) is False
    err = capsys.readouterr()[1]
    assert err == 'E: unable to parse positions in file {}\n'.format(filename)


def test_read_map_raw():
    '''
    Test read_map_raw()
    '''
    indata = ('ali:prot1 OS=Homo sapiens	N	A	C\n' +
              'ali:prot2 OS=Homo sapiens	E	V	A\n' +
              'ali:prot3 OS=Homo sapiens	-	-	-\n' +
              'ali:prot4 OS=Homo sapiens	L	V	P\n' +
              'cons:prot2_human	0.447	0.592	0.3\n' +
              'cons:prot4_human	0.774	0.992	0.4\n' +
              'cons:prot3_human	0.477	0.192	0.5\n' +
              'cons:prot1_human	0.474	0.345	0.6\n' +
              'as:prot1			X\n')
    heads = ['ali:prot1 OS=Homo sapiens',
             'ali:prot2 OS=Homo sapiens',
             'ali:prot3 OS=Homo sapiens',
             'ali:prot4 OS=Homo sapiens']
    seqs = ['NAC', 'EVA', '---', 'LVP']
    data = ['cons:prot2_human	0.447	0.592	0.3',
            'cons:prot4_human	0.774	0.992	0.4',
            'cons:prot3_human	0.477	0.192	0.5',
            'cons:prot1_human	0.474	0.345	0.6',
            'as:prot1			X']
    assert maa.read_map_raw(indata) == (heads, seqs, data)
