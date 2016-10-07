#!/usr/bin/python3

'''
Select accessions from a tabbed blastp result file
The fields are expected to be: "sseqid slen evalue pident length"
'''

import sys

import select_accs_blastp as sab


def test_read_results():
    '''
    Test read_results()
    '''
    rawdata = ('sp|Q11111|PROT_HUMAN	373	0.0	100.000	373\n' +
               'tr|G22222|G22222_PANTR	373	0.0	99.196	373\n' +
               'tr|G33333|G33333_GORGO	373	0.0	98.660	373\n' +
               'tr|G44444|G44444_MACFA	373	0.0	97.855	373\n' +
               'tr|G55555|G55555_MACMU	373	0.0	97.855	373\n' +
               'tr|G66666|G66666_NOMLE	373	0.0	97.855	373\n')

    expected = [['sp|Q11111|PROT_HUMAN', '373', '0.0', '100.000', '373'],
                ['tr|G22222|G22222_PANTR', '373', '0.0', '99.196', '373'],
                ['tr|G33333|G33333_GORGO', '373', '0.0', '98.660', '373'],
                ['tr|G44444|G44444_MACFA', '373', '0.0', '97.855', '373'],
                ['tr|G55555|G55555_MACMU', '373', '0.0', '97.855', '373'],
                ['tr|G66666|G66666_NOMLE', '373', '0.0', '97.855', '373']]

    import tempfile

    filename_blastp = tempfile.mkstemp()[1]
    with open(filename_blastp, 'w') as tmpf:
        tmpf.write(rawdata)

    assert sab.read_results(filename_blastp) == expected


def test_filter_e():
    '''
    Test filter_e()
    '''
    data = [['sp|Q11111|PROT_HUMAN', '373', '0.0', '100.000', '373'],
            ['tr|G22222|G22222_PANTR', '373', '1E-190', '99.196', '373'],
            ['tr|G33333|G33333_GORGO', '373', '1E-150', '78.660', '373'],
            ['tr|G44444|G44444_MACFA', '373', '1E-90', '57.855', '373'],
            ['tr|G55555|G55555_MACMU', '373', '1E-30', '37.855', '373'],
            ['tr|G66666|G66666_NOMLE', '373', '2.0', '17.855', '373']]

    assert sab.filter_e(data) == {'sp|Q11111|PROT_HUMAN'}
    assert sab.filter_e(data, 1E-50) == {'sp|Q11111|PROT_HUMAN',
                                     'tr|G22222|G22222_PANTR',
                                     'tr|G33333|G33333_GORGO',
                                     'tr|G44444|G44444_MACFA'}


def test_main(capsys):
    '''
    Test main()
    '''
    rawdata = ('sp|Q11111|PROT_HUMAN	373	0.0	100.000	373\n' +
               'tr|G22222|G22222_PANTR	373	1E-190	99.196	373\n' +
               'tr|G33333|G33333_GORGO	373	1E-150	78.660	373\n' +
               'tr|G44444|G44444_MACFA	373	1E-90	57.855	373\n' +
               'tr|G55555|G55555_MACMU	373	1E-30	37.855	373\n' +
               'tr|G66666|G66666_NOMLE	373	2.0	17.855	373\n')

    import tempfile

    filename_blastp = tempfile.mkstemp()[1]
    with open(filename_blastp, 'w') as tmpf:
        tmpf.write(rawdata)

    real = 'sp|Q11111|PROT_HUMAN\n'
    sab.main(filename_blastp)
    out, err = capsys.readouterr()
    assert out == real

    expected = {'sp|Q11111|PROT_HUMAN',
                'tr|G22222|G22222_PANTR',
                'tr|G33333|G33333_GORGO',
                ''}
    sab.main(filename_blastp, 1E-140)
    out, err = capsys.readouterr()
    assert set(out.split('\n')) == expected
