#!/usr/bin/env python3
'''
Extract positions from sequences
Reads a FASTA file and a position file
Identifies the sequences, extracts positions,
and prints the merged positions in FASTA format.

Position file:
Identifier\tpos\t...
'''

import sys

import bioinfo

def main(fasta_fname, pos_fname):
    '''
    Reads a FASTA file and a position file
    Identifies the sequences, extracts positions,
    and prints the merged positions in FASTA format.
    '''
    seqdata = bioinfo.read_fasta(fasta_fname)
    posdata = read_posdata(pos_fname)
    for posdat in posdata:
        i = [header[0] for header in enumerate(seqdata[0]) if posdat[0] in header[1]][0]
        print('>' + seqdata[0][i])
        print(''.join([seqdata[1][i][pos-1] for pos in posdat[1:][0]]))


def test_main(capsys):
    '''
    Test main()
    '''
    import tempfile
    data = ('>Protein 1 P12345\n' +
            'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY' +
            'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY' +
            'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n' +
            '>Protein 2 Q12345\n' +
            'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY' +
            'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY' +
            'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n')
    seq_filename = tempfile.mkstemp()[1]
    with open(seq_filename, 'w') as tmpfile:
        tmpfile.write(data)

    data = ('Q12345\t19\t22\t33\t95\n' +
            'P12345\t7\t18\t26\t89\n')

    pos_filename = tempfile.mkstemp()[1]
    with open(pos_filename, 'w') as tmpfile:
        tmpfile.write(data)

    main(seq_filename, pos_filename)
    out = capsys.readouterr()[0]
    expected = ('>Protein 2 Q12345\nWCPR\n' +
                '>Protein 1 P12345\nHVGK\n')
    assert out == expected


def read_posdata(pos_fname):
    '''
    Read positional data
    '''
    with open(pos_fname) as infile:
        data = infile.read().split('\n')
    data = [dat.split('\t') for dat in data if len(dat) > 0]
    data = [[dat[0], [int(pos) for pos in dat[1:]]] for dat in data]
    return data


def test_read_posdata():
    '''
    Test read_posdata()
    '''
    import tempfile
    data = ('P12345\t7\t18\t26\t89\n' +
            'Q12345\t19\t22\t33\t95\n')
    filename = tempfile.mkstemp()[1]
    with open(filename, 'w') as tmpfile:
        tmpfile.write(data)
    expected = [['P12345', [7, 18, 26, 89]],
                ['Q12345', [19, 22, 33, 95]]]

    assert read_posdata(filename) == expected


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {} '.format(sys.argv[0]) +
                         '<FASTA file> <Position file>\n')
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
