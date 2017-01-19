#!/usr/bin/env python3
'''
Filter the sequences, keeping only the ones present in
the provided taxonomy file.
'''

import seq_by_species as sbs


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile
    infasta = ('>sp|Uncharacterized protein OS=Homo sapiens GN=X\n' +
               'ACDEFGHIKLMNPQRSTVWY\n' +
               '>sp|Uncharacterized protein OS=Arabidopsis thaliana GN=X\n' +
               'AAAAAAACCCCCCCCDDDDD\n' +
               '>sp|Uncharacterized protein OS=Felis Cattus GN=X\n' +
               'EEEEEEFFFFFFGGGGGGGG\n')

    filename_fasta = tempfile.mkstemp()[1]
    with open(filename_fasta, 'w') as tmpf:
        tmpf.write(infasta)
    filename_tax = 'tests/testdata/tax_test.tab'
    sbs.main(filename_fasta, filename_tax)
    out = capsys.readouterr()[0]
    expected = ('>sp|Uncharacterized protein OS=Homo sapiens GN=X\n' +
                'ACDEFGHIKLMNPQRSTVWY\n\n' +
                '>sp|Uncharacterized protein OS=Arabidopsis thaliana GN=X\n' +
                'AAAAAAACCCCCCCCDDDDD\n\n')
    assert out == expected


def test_read_species():
    '''
    Test read_species
    '''
    filename_tax = 'tests/testdata/tax_test.tab'
    # quite long list; only sampling
    print(sbs.read_species(filename_tax))
    expected = ['16SrII (Peanut WB group)',
                '16SrXIII (Mexican periwinkle virescence group)',
                '16SrXXII (Nigerian coconut lethal decline (LDN) group)',
                'actinobacterium K18',
                'alpha proteobacterium WTs-6',
                'Bacillus sp. S-(s)-m-D-1(1)',
                'Burkholderiales bacterium RIFCSPLOWO2_12_FULL_61_40',
                'Sclerotinia homoeocarpa mitovirus',
                'Gemmata-like str. Soil9',
                'Naegleria sp. RNG065',
                'Norovirus Hu/GII.4/21859980/AUS/2013',
                'Homo sapiens',
                'Homo sapiens/Rattus norvegicus xenograft',
                'Rattus norvegicus albus',
                'Arabidopsis thaliana x Arabidopsis halleri']
    assert sbs.read_species(filename_tax)[::5] == expected
