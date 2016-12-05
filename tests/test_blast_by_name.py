#!/usr/bin/env python3
'''
Filter the BLAST results based on the last hit with the correct gene name
'''

import blast_by_name as bbn


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    indata = ('>prot1 GN=gene1 SV=1\n' +
              'MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQL\n' +
              'TMVTYGGMAKQPVTASVSMLIFKDLKLRGFWLSQWKKNHSPDEF\n' +
              '>prot2 GN=UNK SV=1\n' +
              'MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVT\n' +
              'FEEFPVALDKALGKLGRQPKQVITF\n' +
              '>prot3 GN=gene1 SV=1\n' +
              'MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVT\n' +
              'FEEFPVALDKALGKLGRQPKQVITF\n'
              '>prot4 GN=gene2 SV=1\n' +
              'MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVT\n' +
              'FEEFPVALDKALGKLGRQPKQVITF\n')
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    expected = ('>prot1 GN=gene1 SV=1\n' +
                'MLVSRRLTGARARAPLLASLLEAWCRQGRTTSSYSAFSEPSHVRALVYGNHGDPAKVIQL\n' +
                'TMVTYGGMAKQPVTASVSMLIFKDLKLRGFWLSQWKKNHSPDEF\n' +
                '>prot2 GN=UNK SV=1\n' +
                'MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVT\n' +
                'FEEFPVALDKALGKLGRQPKQVITF\n' +
                '>prot3 GN=gene1 SV=1\n' +
                'MAALMESVVGRALKFSSTANFRSIRRGETPTLCIKSFSTIMSPPSKAIVYEEHGSPDSVT\n' +
                'FEEFPVALDKALGKLGRQPKQVITF\n')

    bbn.main(filename)
    out = capsys.readouterr()[0]
    assert out == expected
