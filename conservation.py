#!/usr/bin/env python3
'''Calculate the conservation in an alignment'''

import sys

from Bio import AlignIO
from Bio.Align import AlignInfo


def helper_test_getalign():
    '''Generate a FASTA file for use in the testing functions'''
    import tempfile

    indata = '''>gi|0000000|ref|NP_000000.1| Made-up data [Rattus norvegicus]
ACDEFGHIKL
>gi|0000001|ref|NP_000001.1| Made-up data [Arabidopsis thaliana]
ACDEAGHIKL
>gi|0000002|ref|NP_000002.1| Made-up data [Homo Sapiens]
ACDEFGHIKL
>gi|0000003|ref|NP_000003.1| Made-up data [Arabidopsis thaliana]
ADDEEGHILL'''

    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpf:
        tmpf.write(indata)

    return filename


def get_most_conserved(freq_table, align_len):
    '''Determine the most conserved residue and its conservation rate in each position
    Input should be a PSSM object or a list of dicts (PSSM.pssm, but without [0] in each tuple)'''
    result = list()
    # depending on input type
    try:
        num_positions = len(freq_table)
    except TypeError:
        num_positions = len(freq_table.pssm)

    for pos in range(num_positions):
        freq_table_inv = dict((j, i) for i, j in freq_table[pos].items())
        best_conserved = freq_table_inv[max(freq_table_inv)]
        score = freq_table[pos][best_conserved]/align_len
        result.append((score, best_conserved))

    return result


def test_get_most_conserved():
    '''Test get_most_conserved'''
    freq_table = [{'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 4.0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 3.0, 'F': 0, 'D': 1.0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 4.0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 4.0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 1.0,
                   'C': 0, 'F': 2.0, 'D': 0, 'G': 0, 'E': 1.0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 4.0, 'E': 0},
                  {'K': 0, 'H': 4.0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 4.0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 3.0, 'H': 0, 'I': 0, 'L': 1.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0},
                  {'K': 0, 'H': 0, 'I': 0, 'L': 4.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}]

    real = [(1.0, 'A'), (0.75, 'C'), (1.0, 'D'), (1.0, 'E'), (0.5, 'F'),
            (1.0, 'G'), (1.0, 'H'), (1.0, 'I'), (0.75, 'K'), (1.0, 'L')]

    assert get_most_conserved(freq_table, 4) == real


def make_freq_table(alignment):
    '''Make a pssm of an alignment'''
    summary = AlignInfo.SummaryInfo(alignment)
    consensus = summary.dumb_consensus()
    freq_table = summary.pos_specific_score_matrix(consensus)

    return freq_table


def test_make_freq_table():
    '''Test make_freq_table'''
    real = [('A', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 4.0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('C', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 3.0, 'F': 0, 'D': 1.0, 'G': 0, 'E': 0}),
            ('D', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 4.0, 'G': 0, 'E': 0}),
            ('E', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 4.0}),
            ('X', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 1.0,
                   'C': 0, 'F': 2.0, 'D': 0, 'G': 0, 'E': 1.0}),
            ('G', {'K': 0, 'H': 0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 4.0, 'E': 0}),
            ('H', {'K': 0, 'H': 4.0, 'I': 0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('I', {'K': 0, 'H': 0, 'I': 4.0, 'L': 0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('K', {'K': 3.0, 'H': 0, 'I': 0, 'L': 1.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0}),
            ('L', {'K': 0, 'H': 0, 'I': 0, 'L': 4.0, 'A': 0,
                   'C': 0, 'F': 0, 'D': 0, 'G': 0, 'E': 0})]

    alignment = AlignIO.read(helper_test_getalign(), 'fasta')

    assert make_freq_table(alignment).pssm == real


def main(filename, refseq):
    '''Read an alignment in FASTA format
    Calculate the conservation per position'''
    # load alignments
    alignment = AlignIO.read(filename, 'fasta')

    headers = [s.name for s in alignment]
    try:
        refind = headers.index([h for h in headers if refseq in h][0])
    except IndexError:
        sys.stderr.write('E: The reference sequence ({}) not found among the sequences\n'.format(refseq))
        sys.exit(1)

    seqlen = len(str(alignment[refind].seq).replace('-', ''))

    freq_table = make_freq_table(alignment)
    cons = get_most_conserved(freq_table, len(alignment))

    for i in range(len(freq_table.pssm)):
        print(cons[i])


def test_main(capsys):
    '''Test main()'''
    real = ('''(1.0, 'A')\n(0.75, 'C')\n(1.0, 'D')\n(1.0, 'E')\n(0.5, 'F')\n''' +
            '''(1.0, 'G')\n(1.0, 'H')\n(1.0, 'I')\n(0.75, 'K')\n(1.0, 'L')\n''')

    main(helper_test_getalign(), 'NP_000001.1')
    out, err = capsys.readouterr()

    assert out == real


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {0} <alignment file> <reference sequence>\n'.format(sys.argv[0]))
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
