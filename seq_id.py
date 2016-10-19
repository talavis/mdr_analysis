#!/usr/bin/env python3
'''
Calculate the sequence identity between the sequences in an alignment
'''

import sys
import numpy

import bioinfo


def calc_seqids(seqs, skip=True):
    '''
    Calculate the pairwise sequence identities
    Parameter: list of aligned sequences
    Return: list of seq_ids for all comparisons
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
    '''
    # len should be same
    matches = list()
    for i in range(len(sequence1)):
        if skip:
            if sequence1[i] == '-' or sequence2[i] == '-':
                continue
        if sequence1[i] == sequence2[i]:
            matches.append(1)
        else:
            matches.append(0)
    return matches


def main(filename, options=[]):
    '''
    Read a FASTA alignment
    Calculate the pairwise sequence identities
    Optionally print basic statistics
    '''
    config = set_config(options)
    if isinstance(config, str):  # something went wrong
        return config
    heads, seqs = bioinfo.read_fasta(filename)
    seq_ids = calc_seqids(seqs, config[0])
    if config[1]:
        print_stats(seq_ids)


def print_stats(seq_ids):
    '''
    Print the basic stats
    '''
    print('Median: {:.4}'.format(numpy.median(seq_ids)))
    print('Average: {:.4}'.format(numpy.average(seq_ids)))
    print('Max: {:.4}'.format(max(seq_ids)))
    print('Min: {:.4}'.format(min(seq_ids)))


def set_config(options):
    '''
    Change cli options to standardised tuple
    Returns a list with the config if fine
    Returns an error string otherwise
    '''
    accepted = {'gs': (0, True), 'gm': (0, False),
                'sn': (1, False), 'sy': (1, True)}
    config = [True] * (len(accepted)//2)
    for opt in options:
        try:
            tmp = accepted[opt]
            config[tmp[0]] = tmp[1]
        except KeyError:
            return 'Unknown option: {}'.format(opt)
    return config


if __name__ == '__main__':
    if len(sys.argv) < 2:
        USAGE_INFO = ('Usage {0} <alignment file>'.format(sys.argv[0]) +
                      '[gaps: gs/gm stats: sy/sn]\n')
        sys.stderr.write(USAGE_INFO)
        sys.exit()

    RUNSTATUS = main(sys.argv[1])
    if RUNSTATUS is not None:
        sys.stderr.write('E: {}\n'.format(RUNSTATUS))
