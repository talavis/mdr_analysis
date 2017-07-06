#!/usr/bin/env python3
'''
Generate a set of treefiles with sequence names corresponding
to the columns of a FASTA alignment
'''

import sys

import bioinfo

def main(ali_filename, tree_filename):
    '''
    Generate a set of treefiles with sequence names corresponding
    to the columns of a FASTA alignment
    '''
    ali = bioinfo.read_fasta(ali_filename)
    tree_raw = open(tree_filename).read()

    # get rid of path
    startpos = len(ali_filename)-ali_filename[::-1].index('/')
    basename = ali_filename[startpos:ali_filename.index('.', startpos)]
    pos = 0
    for i in range(len(ali[1][0])):
        if ali[1][0][i] == '-':
            continue
        pos += 1
        tree_new = tree_raw
        tree_out = open(basename + str(pos) + '.tree', 'w')
        # assume uniprot headers; accession code between 1st and 2nd |
        # in other words: str[3:str.index('|', 3)
        for head in range(len(ali[0])):
            acc = ali[0][head][3:ali[0][head].index('|', 3)]
            tree_new = tree_new.replace(acc, '{}_{}'.format(ali[1][head][i], acc))
        tree_out.write(tree_new)
        tree_out.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} <FASTA alignment> <tree file>\n'.format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
