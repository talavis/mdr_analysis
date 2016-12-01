#!/usr/bin/env python3
'''
Map ICM positions for active site to conservation data
'''

import map_icm_cons_as as mica


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    # correct
    protdata = ('>sp|Q8N4Q0|ZADH2_HUMAN\n' +
                'MLRLVPTGARAIVDMSYARHFLDFQGSAIPQAMQKLVVTRLSPNFREAVTLSRDCPVPLP\n' +
                'GDGDLLVRNRFVGVNASDINYSAGRYDPSVKPPFDIGFEGIGEVVALGLSASARYTVGQA\n' +
                'VAYMAPGSFAEYTVVPASIATPVPSVKPEYLTLLVSGTTAYISLKELGGLSEGKKVLVTA\n' +
                'AAGGTGQFAMQLSKKAKCHVIGTCSSDEKSAFLKSLGCDRPINYKTEPVGTVLKQEYPEG\n' +
                'VDVVYESVGGAMFDLAVDALATKGRLIVIGFISGYQTPTGLSPVKAGTLPAKLLKKSASV\n' +
                'QGFFLNHYLSKYQAAMSHLLEMCVSGDLVCEVDLGDLSPEGRFTGLESIFRAVNYMYMGK\n' +
                'NTGKIVVELPHSVNSKL\n')
    prot_name = tempfile.mkstemp()[1]
    with open(prot_name, 'w') as tmpf:
        tmpf.write(protdata)

    inasdata = ('- Num  Res. Type ---- SS Molecule ---- Object - sf - sfRatio\n' +
                '   66  asn   Amino    N _  a            2c0c     3.3  0.02  .  a_2c0c.a/^N66\n' +
                '   68  ser   Amino    S H  a            2c0c    17.7  0.14  .  a_2c0c.a/^S68\n' +
                '   71  asn   Amino    N H  a            2c0c    15.8  0.10  .  a_2c0c.a/^N71\n' +
                '   72  tyr   Amino    Y H  a            2c0c    37.0  0.16  .  a_2c0c.a/^Y72\n' +
                '   76  arg   Amino    R _  a            2c0c    47.7  0.20  .  a_2c0c.a/^R76\n' +
                '   77  tyr   Amino    Y _  a            2c0c    96.7  0.41  ~  a_2c0c.a/^Y77\n' +
                '   89  phe   Amino    F _  a            2c0c    28.3  0.13  .  a_2c0c.a/^F89\n' +
                '  106  tyr   Amino    Y _  a            2c0c    20.5  0.09  .  a_2c0c.a/^Y106\n' +
                '  115  met   Amino    M E  a            2c0c    31.3  0.15  .  a_2c0c.a/^M115\n' +
                '  116  ala   Amino    A _  a            2c0c    14.9  0.12  .  a_2c0c.a/^A116\n' +
                '  117  pro   Amino    P _  a            2c0c    85.9  0.57  ~  a_2c0c.a/^P117\n' +
                '  129  ser   Amino    S G  a            2c0c   114.1  0.91  ~  a_2c0c.a/^S129\n' +
                '  130  ile   Amino    I G  a            2c0c    77.1  0.40  ~  a_2c0c.a/^I130\n' +
                '  131  ala   Amino    A G  a            2c0c    11.8  0.10  .  a_2c0c.a/^A131\n' +
                '  132  thr   Amino    T E  a            2c0c    17.6  0.12  .  a_2c0c.a/^T132\n' +
                '  146  val   Amino    V H  a            2c0c    11.8  0.07  .  a_2c0c.a/^V146\n' +
                '  260  ile   Amino    I _  a            2c0c     0.0  0.00  .  a_2c0c.a/^I260\n' +
                '  261  gly   Amino    G _  a            2c0c     0.0  0.00  .  a_2c0c.a/^G261\n' +
                '  262  phe   Amino    F G  a            2c0c    16.9  0.08  .  a_2c0c.a/^F262\n' +
                '  266  tyr   Amino    Y G  a            2c0c     7.5  0.03  .  a_2c0c.a/^Y266\n' +
                '  271  gly   Amino    G _  a            2c0c     0.8  0.01  .  a_2c0c.a/^G271\n' +
                '  272  leu   Amino    L _  a            2c0c    46.7  0.24  .  a_2c0c.a/^L272\n' +
                '  295  phe   Amino    F G  a            2c0c    48.5  0.22  .  a_2c0c.a/^F295\n')
    asdata_name = tempfile.mkstemp()[1]
    with open(asdata_name, 'w') as tmpf:
        tmpf.write(inasdata)

    icmvisdata = ('read pdb "2c0c"\n' +
                  'convertObject a_2c0c. 1==1 no yes yes yes yes yes ' +
                  '""+( 1==2 ? "water=tight ":"" )\n' +
                  'cons = a_2c0c./^N35,^R37,^E38,^A39,^V40,^L42,^C46,^V48,' +
                  '^D55,^L56,^L57,^R59,^N60,^R61,^G64,^I70,^S80,^F89,^G93,' +
                  '^Y106,^Q110,^A111,^S119,^A121,^E122,^V126,^I130,^P135,' +
                  '^K138,^P139,^Y141,^L142,^T143,^L144,^L145,^S147,^G148,' +
                  '^T149,^T150,^A151,^Y152,^I153,^K156,^L158,^S162,^K165,' +
                  '^K166,^A171,^G175,^F179,^M181,^S184,^K186,^C189,^C195,' +
                  '^K200,^A202,^L204,^G208,^C209,^D210,^R211,^P212,^T222,' +
                  '^V223,^Q226,^G231,^V232,^V234,^E237,^R256,^I258,^V259,' +
                  '^S264,^G265,^Q267,^P281,^S288,^S290,^F295,^Y303,^A305,' +
                  '^E312,^V315,^V323,^D324,^G326,^D327\n' +
                  'color a_2c0c. green\n' +
                  'color cons red\n' +
                  'color cool a_2c0c.\n')

    icmvis_name = tempfile.mkstemp()[1]
    with open(icmvis_name, 'w') as tmpf:
        tmpf.write(icmvisdata)

    mica.main(prot_name, '2c0c', icmvis_name, asdata_name)
    out, err = capsys.readouterr()
    print(out, err)
    assert out == ('97\tF\n' +
                   '114\tY\n' +
                   '138\tI\n' +
                   '303\tF\n')


def test_read_vis():
    '''
    Test read_vis()
    '''
    import tempfile

    indata = ('read pdb "2c0c"\n' +
              'convertObject a_2c0c. 1==1 no yes yes yes yes yes ' +
              '""+( 1==2 ? "water=tight ":"" )\n' +
              'cons = a_2c0c./^N35,^R37,^E38,^A39,^V40,^L42,^C46,^V48,' +
              '^D55,^L56,^L57,^R59,^N60,^R61,^G64,^I70,^S80,^F89,^G93\n' +
              'color a_2c0c. green\n' +
              'color cons red\n' +
              'color cool a_2c0c.\n')

    data_name = tempfile.mkstemp()[1]
    with open(data_name, 'w') as tmpf:
        tmpf.write(indata)

    pos = [35, 37, 38, 39, 40, 42, 46, 48, 55, 56,
           57, 59, 60, 61, 64, 70, 80, 89, 93]
    res = ['N', 'R', 'E', 'A', 'V', 'L', 'C', 'V', 'D', 'L',
           'L', 'R', 'N', 'R', 'G', 'I', 'S', 'F', 'G']
    assert mica.read_vis(data_name) == (pos, res)
