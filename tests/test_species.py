#!/usr/bin/env python3
'''
Tests for the species module
'''

import species
import bioinfo
import taxnode


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
    '''Test make_taxtree()'''
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
    '''Test match_taxonomy()'''
    species_tup = ('Rattus norvegicus', 'Arabidopsis thaliana', 'Homo sapiens')
    filename = 'testdata/tax_test.tab'
    assert [col[0] for col in species.match_taxonomy(filename, species_tup)] == ['9606', '10116', '3702']


def test_print_tax_tree(capsys):
    '''Test print_tax_tree()'''
    specs = [[1]*8+['a; a1; a11'],
             [1]*8+['a; a1; a12'],
             [1]*8+['a; a2; a21'],
             [1]*8+['b; b1; b11']]

    root = species.make_taxtree(specs)
    species.print_tax_tree(root)
    out, err = capsys.readouterr()
    expected = ('\\root\n  \\a\n    \\a1\n      \\a11\n      \\a12\n    ' +
                '\\a2\n      \\a21\n  \\b\n    \\b1\n      \\b11\n')
    assert out == expected


def test_main(capsys):
    '''Test main()'''
    import tempfile
    infasta = ('>sp|Uncharacterized protein OS=Homo sapiens GN=X\n' +
               '>sp|Uncharacterized protein OS=Arabidopsis thaliana GN=X\n' +
               '>sp|Uncharacterized protein OS=Rattus norvegicus GN=X\n')
    filename_fasta = tempfile.mkstemp()[1]
    with open(filename_fasta, 'w') as tmpf:
        tmpf.write(infasta)
    filename_tax = 'testdata/tax_test.tab'
    species.main(filename_fasta, filename_tax)
    out, err = capsys.readouterr()
    expected = ('\\root\n  \\Eukaryota\n    \\Metazoa\n      \\Chordata\n        ' +
                '\\Craniata\n          \\Vertebrata\n            \\Euteleostomi\n' +
                '              \\Mammalia\n                \\Eutheria\n' +
                '                  \\Euarchontoglires\n                    ' +
                '\\Primates\n                      \\Haplorrhini\n' +
                '                        \\Catarrhini\n                          ' +
                '\\Hominidae\n                            \\Homo\n                    ' +
                '\\Glires\n                      \\Rodentia\n                        ' +
                '\\Sciurognathi\n                          \\Muroidea\n' +
                '                            \\Muridae\n                              ' +
                '\\Murinae\n                                \\Rattus\n    ' +
                '\\Viridiplantae\n      \\Streptophyta\n        \\Embryophyta\n' +
                '          \\Tracheophyta\n            \\Spermatophyta\n              ' +
                '\\Magnoliophyta\n                \\eudicotyledons\n                  ' +
                '\\Gunneridae\n                    \\Pentapetalae\n' +
                '                      \\rosids\n                        \\malvids\n' +
                '                          \\Brassicales\n                            ' +
                '\\Brassicaceae\n                              \\Camelineae\n' +
                '                                \\Arabidopsis\n' +
                'Highest common: Eukaryota\n')

    assert out == expected
