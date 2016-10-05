#!/usr/bin/python3

'''
Select accessions from a tabbed blastp result file
The fields are expected to be: "sseqid slen evalue pident length"
'''

import sys


def read_results(filename):
    '''
    Read a blastp file
    Returns a list of hits
    '''

    data = list()

    with open(filename) as infile:
        for line in infile:
            data.append(line.rstrip().split('\t'))

    return data


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

    assert read_results(filename_blastp) == expected


def filter_e(results, evalue=0.0):
    '''
    Filter the results based on E values
    Returns a set of accessions (column 0)
    '''
    matches = set()
    for res in results:
        if float(res[2]) <= evalue:
            matches.add(res[0])
    return matches


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

    assert filter_e(data) == {'sp|Q11111|PROT_HUMAN'}
    assert filter_e(data, 1E-50) == {'sp|Q11111|PROT_HUMAN',
                                     'tr|G22222|G22222_PANTR',
                                     'tr|G33333|G33333_GORGO',
                                     'tr|G44444|G44444_MACFA'}


def main(filename, evalue=0.0):
    '''
    Read a blastp file and filter the results based on E values
    Prints the resulting set of accessions
    '''
    data = read_results(filename)
    accs = filter_e(data, evalue)

    for acc in accs:
        print(acc)


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
    main(filename_blastp)
    out, err = capsys.readouterr()
    assert out == real

    expected = {'sp|Q11111|PROT_HUMAN',
                'tr|G22222|G22222_PANTR',
                'tr|G33333|G33333_GORGO',
                ''}
    main(filename_blastp, 1E-140)
    out, err = capsys.readouterr()
    assert set(out.split('\n')) == expected


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.stderr.write('Usage: {} '.format(sys.argv[0]) +
                         '<blastp file> ' +
                         '[e value limit]\n')
        sys.exit(1)
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main(sys.argv[1], float(sys.argv[2]))
