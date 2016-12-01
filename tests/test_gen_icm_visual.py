#!/usr/bin/env python3
'''
Generate the commands needed for visualisation of data on a structure in ICM
'''

import gen_icm_visual as giv


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile

    # correct
    inseq = ('>prot1\n' +
             'MSTAGKVIKCKAAVLW\n')
    fasta_name = tempfile.mkstemp()[1]
    with open(fasta_name, 'w') as tmpf:
        tmpf.write(inseq)

    inpos = '2\n4\n7\n'
    pos_name = tempfile.mkstemp()[1]
    with open(pos_name, 'w') as tmpf:
        tmpf.write(inpos)

    expected = ('read pdb "1u3w"\n'
                'convertObject a_1u3w. 1==1 no yes yes yes yes yes ' +
                '""+( 1==2 ? "water=tight ":"" )\n' +
                'cons = a_1u3w./^S1,^A3,^V6\n' +
                'color a_1u3w. green\n' +
                'color cons red\n' +
                'cool a_1u3w.\n')
    giv.main(fasta_name, '1u3w', pos_name)
    out = capsys.readouterr()[0]
    assert out == expected

    # incorrect protein seq file name
    assert giv.main('__incorrect', '1u3w', pos_name) is False
    err = capsys.readouterr()[1]
    expected = 'E: Cannot find FASTA file: __incorrect\n'
    assert err == expected

    # no protein sequence
    inseq = ('>prot1\n' +
             '\n')
    fasta_name = tempfile.mkstemp()[1]
    with open(fasta_name, 'w') as tmpf:
        tmpf.write(inseq)
    assert giv.main('__incorrect', '1u3w', pos_name) is False
    capsys.readouterr()

    # incorrect structure
    inseq = ('>prot1\n' +
             'MSTAGKVIKCKAAVLW\n')
    fasta_name = tempfile.mkstemp()[1]
    with open(fasta_name, 'w') as tmpf:
        tmpf.write(inseq)
    assert giv.main(fasta_name, '__asdfgh', pos_name) is False
    err = capsys.readouterr()[1]
    expected = 'E: could not retrieve sequence for structure __asdfgh\n'
    assert err == expected

    # protein and structure seq not matching
    inseq = ('>prot1\n' +
             'MSTAGKVIKCKAAVLW\n')
    fasta_name = tempfile.mkstemp()[1]
    with open(fasta_name, 'w') as tmpf:
        tmpf.write(inseq)
    assert giv.main(fasta_name, '4ins', pos_name) is False
    err = capsys.readouterr()[1]
    expected = 'E: the sequence STAGKV is not found in the second protein\n'
    assert err == expected

    # incorrect positions
    inpos = '500\n'
    pos_name = tempfile.mkstemp()[1]
    with open(pos_name, 'w') as tmpf:
        tmpf.write(inpos)
    assert giv.main(fasta_name, '1u3w', pos_name) is False
    err = capsys.readouterr()[1]
    expected = ('E: requested positions outside the protein sequence\n' +
                'E: position outside protein (499)\n')
    assert err == expected


def test_make_icm_res():
    '''
    Test make_icm_res()
    '''
    assert giv.make_icm_res('A', 3) == '^A4'
