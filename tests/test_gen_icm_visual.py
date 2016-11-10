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

    expected = ('cons = {' +
                'a_1u3w./^S1 ' +
                'a_1u3w./^A3 ' +
                'a_1u3w./^V6' +
                '}\n')
    giv.main(fasta_name, '1u3w', pos_name)
    out = capsys.readouterr()[0]
    assert out == expected
    

def test_make_icm_res():
    '''
    Test make_icm_res()
    '''
    assert giv.make_icm_res('A', 3, '1u3w') == 'a_1u3w./^A4'
