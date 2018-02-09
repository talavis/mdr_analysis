#!/usr/bin/env python3
'''
Calculate the sequence identity between the sequences in an alignment
'''

import sys
from collections import namedtuple

import numpy

import bioinfo


class DifferentLengthsError(Exception):
    '''
    Used when two sequences have different lengths
    '''
    pass


class IncorrectOptionError(Exception):
    '''
    Used when an incorrect option is provided
    '''
    pass


def calc_seqids(seqs, skip=True):
    '''
    Calculate the pairwise sequence identities of a group of sequences

    Args:
        seqs (list, tuple): aligned sequences as strings
        skip (bool) = skip comparison of positions with gaps

    Returns:
        list: sequence identities for all comparisons
    '''
    seq_ids = list()
    for i in range(len(seqs)-1):
        for j in range(i+1, len(seqs)):
            pos_id = compare_res(seqs[i], seqs[j], skip)
            seq_id = sum(pos_id)/len(pos_id)
            seq_ids.append(seq_id)
    return seq_ids


def compare_res(sequence1, sequence2, skip=True):
    '''
    Compare the residues in two sequences

    Args:
        sequence1 (str): the first sequence to compare
        sequence2 (str): the second sequence to compare
        skip (bool) = skip comparison of positions with gaps

    Returns:
        list: position match information
    '''
    if len(sequence1) != len(sequence2):
        raise DifferentLengthsError('The sequences have different lengths')

    matches = list()
    for pos in zip(sequence1, sequence2):
        if skip:
            if pos[0] == '-' or pos[1] == '-':
                continue
        if pos[0] == pos[1]:
            if pos[0] == '-' and pos[1] == '-':
                continue
            matches.append(1)
        else:
            matches.append(0)
    return matches


def main(filename, options=()):
    '''
    Read a FASTA alignment
    Calculate the pairwise sequence identities
    Optionally print basic statistics

    Args:
        filename (str): filename of the FASTA alignment file
        options (list): options, see set_config
    '''
    config = set_config(options)

    _, seqs = bioinfo.read_fasta(filename)
    seq_ids = calc_seqids(seqs, config.skip_gaps)
    if config.print_stats:
        print_stats(seq_ids)


def print_stats(seq_ids):
    '''
    Print the basic stats of a group of comparisons

    Args:
        seq_ids (list, tuple): sequence identities
    '''
    print('Median: {:.4}'.format(numpy.median(seq_ids)))
    print('Average: {:.4}'.format(numpy.average(seq_ids)))
    print('Max: {:.4}'.format(max(seq_ids)))
    print('Min: {:.4}'.format(min(seq_ids)))


def set_config(options):
    '''
    Generate a configuration based on the given options

    Args:
        options (list, tuple): configuration options

    Returns:
        namedtuple: configuration
    '''
    Config = namedtuple('Config', ['skip_gaps', 'print_stats'])
    accepted = {'gm': (0, False), 'gs': (0, True),
                'sn': (1, False), 'sy': (1, True)}
    tmp_conf = [True]*(len(accepted)//2)
    for opt in options:
        if opt in accepted:
            tmp = accepted[opt]
            tmp_conf[tmp[0]] = tmp[1]
        else:
            raise IncorrectOptionError(f'Unknown option: {opt}')
    config = Config(skip_gaps=tmp_conf[0], print_stats=tmp_conf[1])
    return config


if __name__ == '__main__':
    if len(sys.argv) < 2:
        USAGE_INFO = ('Usage {0} <alignment file>'.format(sys.argv[0]) +
                      '[gaps: gs/gm stats: sy/sn]\n')
        sys.stderr.write(USAGE_INFO)
        sys.exit()

    main(sys.argv[1], sys.argv[2:])
