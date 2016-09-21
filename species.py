#!/usr/bin/env python3
'''Calculating species memberships for sequences in a FASTA file'''

import sys

import bioinfo
import taxnode


def find_highest_common(node):
    '''
    Find the first branching of the tree, counting from the given node
    Return: name of the branching node
    If there are no branches, return the name of the node at the end of the tree
    '''
    while len(node.children) == 1:
        node = node.children[0]
    return node.name


def test_find_highest_common():
    '''
    Test find_highest_common()
    '''
    specs1 = [[1]*8+['a; a1; a11'],
              [1]*8+['a; a1; a12'],
              [1]*8+['a; a2; a21'],
              [1]*8+['b; b1; b11']]
    root = make_taxtree(specs1)
    assert find_lowest_common(root) == 'root'
    specs2 = [[1]*8+['a; a1; a11'],
              [1]*8+['a; a1; a12'],
              [1]*8+['a; a1; a13'],
              [1]*8+['a; a1; b14']]
    root = make_taxtree(specs2)
    assert find_highest_common(root) == 'a1'


def make_taxtree(hits):
    '''Generate a taxonomy tree from a taxonomy list
    returns the root node'''
    root = taxnode.TaxNode('root')
    for hit in range(len(hits)):
        taxonomy = [species.strip() for species in hits[hit][8].split(';')]
        node = root
        for taxa in range(len(taxonomy)):
            if taxonomy[taxa] not in [child.name for child in node.children]:
                node = taxnode.TaxNode(taxonomy[taxa], node)
            else:
                node = [child for child in node.children if child.name == taxonomy[taxa]][0]
    return root


def test_make_taxtree():
    '''Test make_taxtree()'''
    specs = [[1]*8+['a; a1; a11'],
             [1]*8+['a; a1; a12'],
             [1]*8+['a; a2; a21'],
             [1]*8+['b; b1; b11']]

    root = make_taxtree(specs)
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


def match_taxonomy(filename, species):
    '''Match a species list with the taxonomy
    filename: filename for a UniProt taxonomy file'''
    with open(filename) as infile:
        lines = [line for line in infile.read().split('\n') if len(line) > 0]
        hits = [line.split('\t') for line in lines if line.split('\t')[2] in species]

    return hits


def test_match_taxonomy():
    '''Test match_taxonomy()'''
    species = ('Rattus norvegicus', 'Arabidopsis thaliana', 'Homo sapiens')
    filename = 'tests/tax_test.tab'
    assert [col[0] for col in match_taxonomy(filename, species)] == ['9606', '10116', '3702']


def print_tax_tree(node, level=0):
    '''Print a taxonomy tree, starting from the root'''
    print('{}\\{}'.format('  '*level, node.name))
    children = tuple(node.children)
    if len(children) > 0:
        for i in range(len(children)):
            print_tax_tree(children[i], level+1)


def test_print_tax_tree(capsys):
    '''Test print_tax_tree()'''
    specs = [[1]*8+['a; a1; a11'],
             [1]*8+['a; a1; a12'],
             [1]*8+['a; a2; a21'],
             [1]*8+['b; b1; b11']]

    root = make_taxtree(specs)
    print_tax_tree(root)
    out, err = capsys.readouterr()
    expected = ('\\root\n  \\a\n    \\a1\n      \\a11\n      \\a12\n    ' +
                '\\a2\n      \\a21\n  \\b\n    \\b1\n      \\b11\n')
    assert out == expected


def main(fasta_name, taxonomy_name):
    '''Calculating species memberships for sequences in a FASTA file'''
    headers, sequences = bioinfo.read_fasta(fasta_name)
    species = [bioinfo.get_species(h) for h in headers]

    tax_hits = match_taxonomy(taxonomy_name, species)
    tax_tree = make_taxtree(tax_hits)
    print_tax_tree(tax_tree)
    print('Highest common: {}'.format(find_highest_common(tax_tree)))


def test_main(capsys):
    '''Test main()'''
    import tempfile
    infasta = ('>sp|Uncharacterized protein OS=Homo sapiens GN=X\n' +
               '>sp|Uncharacterized protein OS=Arabidopsis thaliana GN=X\n' +
               '>sp|Uncharacterized protein OS=Rattus norvegicus GN=X\n')
    filename_fasta = tempfile.mkstemp()[1]
    with open(filename_fasta, 'w') as tmpf:
        tmpf.write(infasta)
    filename_tax = 'tests/tax_test.tab'
    main(filename_fasta, filename_tax)
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


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <FASTA file> <taxonomy file>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
