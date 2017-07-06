#!/usr/bin/env python3
'''
Calculating species memberships for sequences in a FASTA file
'''

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


def make_taxtree(hits):
    '''
    Generate a taxonomy tree from a taxonomy list
    returns the root node
    '''
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


def match_taxonomy(filename, species):
    '''
    Match a species list with the taxonomy
    filename: filename for a UniProt taxonomy file
    '''
    with open(filename) as infile:
        lines = [line for line in infile.read().split('\n') if line]
        hits = [line.split('\t') for line in lines if line.split('\t')[2] in species]

    return hits


def print_tax_tree(node, level=0):
    '''
    Print a taxonomy tree, starting from the root
    '''
    print('{}\\{}\t{}'.format('  '*level, node.name, node.child_end_count()))
    children = tuple(node.children)
    if children:
        for i in range(len(children)):
            print_tax_tree(children[i], level+1)


def main(fasta_name, taxonomy_name, options=None):
    '''
    Calculating species memberships for sequences in a FASTA file
    '''
    headers, sequences = bioinfo.read_fasta(fasta_name)
    species = [bioinfo.get_species(h) for h in headers]

    tax_hits = match_taxonomy(taxonomy_name, species)
    tax_tree = make_taxtree(tax_hits)

    print_tree = True
    if options:
        if options[0] == 'or':
            print_tree = False
    if print_tree:
        print_tax_tree(tax_tree)
    print('Highest common: {}'.format(find_highest_common(tax_tree)))


if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.stderr.write('Usage: {} <FASTA file> '.format(sys.argv[0]) +
                         '<taxonomy file> [only print root: or]\n')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
