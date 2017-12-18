#!/usr/bin/env python3
'''
Calculate the conservation in an alignment
'''

import sys

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped


def conservation(filename, refseq=None, group_res=False, ign_gaps=False):
    '''
    Read an alignment in FASTA format
    Calculate the conservation per position
    Input: filename - filename of alignment in FASA format
    refseq - only keep positions where sequence refseq has a residue
    group_res - group the residues by residue type
    ign_gaps - ignore gaps when calculating conservation
    '''
    alignment = AlignIO.read(filename, 'fasta')

    headers = [s.name for s in alignment]
    if refseq:
        try:
            refind = headers.index([h for h in headers if refseq in h][0])
        except IndexError:
            error = ('E: The reference sequence ({}) ' +
                     'not found among the sequences\n').format(refseq)
            sys.stderr.write(error)
            return False

    if group_res:
        sequences = [str(ali.seq) for ali in alignment]
        sequences = group_res_prop(sequences)
        for i in range(len(alignment)):
            alignment[i].seq = Seq(sequences[i], Gapped(IUPAC.protein, '-'))

    freq_table = make_freq_table(alignment)
    cons = get_most_conserved(freq_table, ign_gaps)

    if refseq is None:
        refseq_p = ' '*len(alignment[0].seq)
    else:
        refseq_p = str(alignment[refind].seq)
        print('# {}'.format(refseq))
    pos = 1
    for i in range(len(freq_table.pssm)):
        if refseq_p[i] != '-':
            print('{ps}\t{rs}\t{mc}\t{rate:.3}'.format(ps=pos,
                                                       rs=refseq_p[i],
                                                       mc=cons[i][1],
                                                       rate=cons[i][0]))
            pos += 1


def get_most_conserved(freq_table, ign_gaps=False):
    '''
    Determine the most conserved residue and its conservation
    rate in each position.
    ign_gaps - ignore gaps in the calculation of the conservation
    Input: a PSSM frequency table or a list of dicts
    (PSSM.pssm, but without [0] in each tuple)
    Return: a list containing (score, residue)
    '''
    # depending on input type
    try:
        # dictionary with counts
        num_positions = len(freq_table)
    except TypeError:
        # frequency table from Biopython
        num_positions = len(freq_table.pssm)
    result = [0] * num_positions
    for pos in range(num_positions):
        if ign_gaps:
            freq_table[pos]['-'] = 0
        align_size = sum(freq_table[pos].values())
        values = list(freq_table[pos].values())
        if values.count(max(values)) == 1:
            freq_table_inv = dict((j, i) for i, j in freq_table[pos].items())
            best_conserved = freq_table_inv[max(freq_table_inv)]
            score = freq_table[pos][best_conserved]/align_size
        else:
            best_conserved = 'X'
            try:
                score = max(values)/align_size
            except ZeroDivisionError:
                # only gaps in positions; shouldn't occur in normal alignments
                score = 0.0
        result[pos] = (score, best_conserved)

    return result


def group_res_prop(sequences):
    '''
    Transform an alignment by grouping residues by their properties:
    I - ILV
    D - DE
    K - KR
    S - ST
    N - NQ
    Input: sequences - list of protein sequences
    Return: the same alignment with the relevant residues replaced
    by their group names
    '''
    for i in range(len(sequences)):
        sequences[i] = sequences[i].replace('L', 'I')
        sequences[i] = sequences[i].replace('V', 'I')
        sequences[i] = sequences[i].replace('E', 'D')
        sequences[i] = sequences[i].replace('R', 'K')
        sequences[i] = sequences[i].replace('Q', 'N')
        sequences[i] = sequences[i].replace('T', 'S')
    return sequences


def make_freq_table(alignment):
    '''
    Make a pssm of an alignment
    Input: BioPython alignment object
    Return: BioPython frequency list [(cons_res, {res:freq})]
    '''
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus()
    freq_table = summary.pos_specific_score_matrix(consensus)

    return freq_table


def parse_parameters(params):
    '''
    Parse the commandline parameters
    '''
    filename = params[0]
    group_res = False
    reference = None
    ign_gaps = False
    for param in params[1:]:
        if 'group_res=' in param:
            choice = param[param.index('=')+1:]
            if choice.lower() == 'y':
                group_res = True
            elif choice.lower() != 'n':
                raise ValueError('E: incorrect parameter ({})'.format(param))
        elif 'reference=' in param:
            reference = param[param.index('=')+1:]
        elif 'ign_gaps=' in param:
            choice = param[param.index('=')+1:]
            if choice.lower() == 'y':
                ign_gaps = True
            elif choice.lower() != 'n':
                raise ValueError('E: incorrect parameter ({})'.format(param))
        else:
            raise ValueError('E: incorrect parameter ({})'.format(param))
    return (filename, reference, group_res, ign_gaps)


def print_use(base):
    '''
    Print the command line usage
    '''
    sys.stderr.write(('Usage: {0} '.format(base) +
                      '<alignment file> ' +
                      '[reference=seq] ' +
                      '[group_res=y/N] ' +
                      '[ign_gaps=y/N]\n'))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print_use(sys.argv[0])
        sys.exit(1)
    try:
        conservation(*parse_parameters(sys.argv[1:]))
    except ValueError as err:
        print(err)
        sys.exit(1)
