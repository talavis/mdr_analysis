#!/usr/bin/env python3
'''Calculating species memberships for sequences in a FASTA file'''

import sys

import bioinfo
import taxnode

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

# def test_make_taxtree():
#     '''Test make_taxtree()'''
#     pass

def match_taxonomy(filename, species):
    '''Match a species list with the taxonomy
    filename: filename for a UniProt taxonomy file'''
    with open(filename) as infile:
        lines = [line for line in infile.read().split('\n') if len(line) > 0]
        hits = [line.split('\t')  for line in lines if line.split('\t')[2] in species]

    return hits

# def test_match_taxonomy():
#     '''Test match_taxonomy()'''
#     pass

def print_tax_tree(node, level=0):
    '''Print a taxonomy tree, starting from the root'''
    print('{}\\{}'.format(' '*level, node.name))
    children = tuple(node.children)
    if len(children) > 0:
        for i in range(len(children)):
            print_tax_tree(children[i], level+1)

# def test_print_tax_tree():
#     '''Test print_tax_tree()'''
#     pass

def main(fasta_name, taxonomy_name):
    '''Calculating species memberships for sequences in a FASTA file'''
    headers, sequences = bioinfo.read_fasta(fasta_name)
    species = [bioinfo.get_species(h) for h in headers]

    tax_hits = match_taxonomy(taxonomy_name, species)
    tax_tree = make_taxtree(tax_hits)
    print_tax_tree(tax_tree)

# def test_main():
#     '''Test main()'''
#     pass

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <FASTA file> <taxonomy file>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
