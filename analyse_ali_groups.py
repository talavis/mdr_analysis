#!/usr/bin/env python3
'''
Analyse a protein sequence alignment based on properties rather
than exact residues. 
'''

import sys

import bioinfo

def analyse_ali(sequences, groups = True):
    '''
    Analyse an alignment based on conservation.
    Apart from only analysing exact residues,
    also allow analysis based on property groups,
    e.g. hydrophobic.
    '''


def test_analyse_ali():
    '''
    Test analyse_ali()
    '''


def main(filename):
    '''
    Read a FASTA alignment and perform a conservation analysis
    using groups of residues
    '''


def test_main(capsys):
    '''
    Test main()
    '''


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: {} <FASTA alignment>\n'format(sys.argv[0]))
        sys.exit(1)
    main(sys.argv[1])
