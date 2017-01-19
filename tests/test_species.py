#!/usr/bin/env python3
'''
Tests for the species module
'''

import species


def test_find_highest_common():
    '''
    Test find_highest_common()
    '''
    specs1 = [[1]*8+['a; a1; a11'],
              [1]*8+['a; a1; a12'],
              [1]*8+['a; a2; a21'],
              [1]*8+['b; b1; b11']]
    root = species.make_taxtree(specs1)
    assert species.find_highest_common(root) == 'root'
    specs2 = [[1]*8+['a; a1; a11'],
              [1]*8+['a; a1; a12'],
              [1]*8+['a; a1; a13'],
              [1]*8+['a; a1; b14']]
    root = species.make_taxtree(specs2)
    assert species.find_highest_common(root) == 'a1'


def test_make_taxtree():
    '''
    Test make_taxtree()
    '''
    specs = [[1]*8+['a; a1; a11'],
             [1]*8+['a; a1; a12'],
             [1]*8+['a; a2; a21'],
             [1]*8+['b; b1; b11']]

    root = species.make_taxtree(specs)
    assert root.parent is None
    assert set(child.name for child in root.children) == {'a', 'b'}
    node = [child for child in root.children if child.name == 'a'][0]
    assert node.parent is root
    assert set(child.name for child in node.children) == {'a1', 'a2'}
    node2 = [child for child in node.children if child.name == 'a1'][0]
    assert node2.parent is node
    assert set(child.name for child in node2.children) == {'a11', 'a12'}
    node2 = [child for child in node.children if child.name == 'a2'][0]
    assert node2.parent is node
    assert set(child.name for child in node2.children) == {'a21'}
    node = [child for child in root.children if child.name == 'b'][0]
    assert node.parent is root
    assert set(child.name for child in node.children) == {'b1'}
    node2 = [child for child in node.children if child.name == 'b1'][0]
    assert node2.parent is node
    assert set(child.name for child in node2.children) == {'b11'}


def test_match_taxonomy():
    '''
    Test match_taxonomy()
    '''
    species_tup = ('Rattus norvegicus', 'Arabidopsis thaliana', 'Homo sapiens')
    filename = 'tests/testdata/tax_test.tab'
    result = species.match_taxonomy(filename, species_tup)
    ids = ['9606', '10116', '3702']
    assert [col[0] for col in result] == ids


def test_print_tax_tree(capsys):
    '''
    Test print_tax_tree()
    '''
    specs = [[1]*8+['a; a1; a11'],
             [1]*8+['a; a1; a12'],
             [1]*8+['a; a2; a21'],
             [1]*8+['b; b1; b11']]

    root = species.make_taxtree(specs)
    species.print_tax_tree(root)
    out = capsys.readouterr()[0]
    expected = ('\\root\t4\n  \\a\t3\n    \\a1\t2\n      \\a11\t1\n      ' +
                '\\a12\t1\n    \\a2\t1\n      \\a21\t1\n  \\b\t1\n    ' +
                '\\b1\t1\n      \\b11\t1\n')
    assert out == expected


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile
    infasta = ('>sp|Uncharacterized protein OS=Homo sapiens GN=X\n' +
               '>sp|Uncharacterized protein OS=Arabidopsis thaliana GN=X\n' +
               '>sp|Uncharacterized protein OS=Rattus norvegicus GN=X\n')
    filename_fasta = tempfile.mkstemp()[1]
    with open(filename_fasta, 'w') as tmpf:
        tmpf.write(infasta)
    filename_tax = 'tests/testdata/tax_test.tab'
    species.main(filename_fasta, filename_tax)
    out = capsys.readouterr()[0]
    expected = ('\\root\t3\n  \\Eukaryota\t3\n    \\Metazoa\t2\n      ' +
                '\\Chordata\t2\n        \\Craniata\t2\n          ' +
                '\\Vertebrata\t2\n            \\Euteleostomi\t2\n' +
                '              \\Mammalia\t2\n                \\Eutheria' +
                '\t2\n                  \\Euarchontoglires\t2\n         ' +
                '           \\Primates\t1\n                      ' +
                '\\Haplorrhini\t1\n                        \\Catarrhini\t1\n' +
                '                          \\Hominidae\t1\n' +
                '                            \\Homo\t1\n    ' +
                '                \\Glires\t1\n                      ' +
                '\\Rodentia\t1\n                        \\Sciurognathi\t1\n' +
                '                          \\Muroidea\t1\n                 ' +
                '           \\Muridae\t1\n                              ' +
                '\\Murinae\t1\n                                \\Rattus\t1\n' +
                '    \\Viridiplantae\t1\n      \\Streptophyta\t1\n        ' +
                '\\Embryophyta\t1\n          \\Tracheophyta\t1\n            ' +
                '\\Spermatophyta\t1\n              \\Magnoliophyta\t1\n' +
                '                \\eudicotyledons\t1\n                  ' +
                '\\Gunneridae\t1\n                    \\Pentapetalae\t1\n' +
                '                      \\rosids\t1\n                        ' +
                '\\malvids\t1\n                          \\Brassicales\t1\n' +
                '                            \\Brassicaceae\t1\n' +
                '                              \\Camelineae\t1\n' +
                '                                \\Arabidopsis\t1\n' +
                'Highest common: Eukaryota\n')
    assert out == expected

    # test parameter "or"
    filename_fasta = tempfile.mkstemp()[1]
    with open(filename_fasta, 'w') as tmpf:
        tmpf.write(infasta)
    filename_tax = 'tests/testdata/tax_test.tab'
    species.main(filename_fasta, filename_tax, ['or'])
    out = capsys.readouterr()[0]
    assert out == 'Highest common: Eukaryota\n'
