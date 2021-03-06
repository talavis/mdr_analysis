#!/usr/bin/env python3
'''
Map active site data from ICM to an alignment or an alignment data map
Use the output of Res(site_def) from ICM
'''

import sys

import bioinfo


def find_refprot_index(refprot, headers):
    '''
    Find the index of the reference sequence
    '''
    try:
        refhead = [head for head in headers if refprot in head]
        refind = headers.index(refhead[0])
    except IndexError:
        error = 'E: reference sequence {} not found\n'.format(refprot)
        sys.stderr.write(error)
        return False

    return refind


def main(map_name, refprot, struct_name, icmres_file):
    '''
    Map active site data to an alignment
    '''
    inmap = open(map_name).read()
    if inmap[0] == '>':
        # FASTA
        heads, seqs = bioinfo.read_fasta_raw(inmap)
        data = []
    else:
        # data map
        heads, seqs, data = read_map_raw(inmap)

    structseq = bioinfo.get_structseq(struct_name)
    if structseq is False:
        return False
    structseq = structseq[1][0]

    as_data = read_icmdata(icmres_file)
    if as_data is False:
        return False
    # shift -1 on positions
    as_data[0] = [pos-1 for pos in as_data[0]]
    if as_data is False:
        return False
    refind = find_refprot_index(refprot, heads)
    if refind is False:
        return False
    positions = map_sequences(structseq, seqs[refind], as_data[0], as_data[1])
    if positions is False:
        return False

    asali_map = map_as(heads, seqs, refprot, positions)

    for i in range(len(seqs)):
        if heads[0][:4] == 'ali:':
            print('{}\t{}'.format(heads[i], '\t'.join([c for c in seqs[i]])))
        else:
            print('ali:{}\t{}'.format(heads[i], '\t'.join([c for c in seqs[i]])))
    for i in range(len(data)):
        print(data[i])
    print('{}\t{}'.format('as:{}_{}'.format(refprot, struct_name),
                          '\t'.join(asali_map)))


def map_as(headers, seqs, refprot, as_pos):
    '''
    Map active site positions to an alignment
    '''
    refind = find_refprot_index(refprot, headers)
    if refind is False:
        return False
    ali_data = ['']*len(seqs[0])
    for i in range(len(as_pos)):
        pos = map_pos(as_pos[i]+1, seqs[refind])
        ali_data[pos] = 'X'
    return ali_data


def map_pos(position, sequence):
    '''
    Map a position in a protein to the position in the alignment
    input: pos1 = 1
    output: pos1 = 0
    '''
    pos = 0
    for i in range(len(sequence)):
        if sequence[i] != '-':
            pos += 1
        if pos == position:
            return i
    err = 'E: position {} not found in {}\n'.format(position, sequence[:20])
    sys.stderr.write(err)
    return False


def map_sequences(query_seq, ref_seq, query_pos, query_res):
    '''
    Map the positions in one sequence to those in another sequence
    Intended for structure vs protein sequence and vice versa
    pos 1 = 0
    '''
    ref_seq = ref_seq.replace('-', '')
    new_pos = [0] * len(query_pos)
    for pos in range(len(query_pos)):
        ind = query_pos[pos]
        if ind >= len(query_seq):
            error = 'E: position outside protein ({})\n'.format(ind)
            sys.stderr.write(error)
            return False
        if query_seq[ind] != query_res[pos]:
            error = ('E: the protein sequence does not match the position data; ' +
                     'Position {} should be {}, but is {}\n'.format(query_pos[pos]+1,
                                                                    query_res[pos],
                                                                    query_seq[ind]))
            sys.stderr.write(error)
        try:
            new_pos[pos] = ref_seq.index(query_seq[ind:ind+6])
        except ValueError:
            error = ('E: the sequence {} '.format(query_seq[ind:ind+6]) +
                     'is not found in the second protein\n')
            sys.stderr.write(error)
            return False
    return new_pos


def read_icmdata(filename):
    '''
    Read active site data
    Should be the output of Res(site_def) from ICM
    '''
    try:
        with open(filename) as infile:
            raw = infile.read()
        lines = raw.split('\n')
    except FileNotFoundError:
        error = 'E: file {} not found\n'.format(filename)
        sys.stderr.write(error)
        return False
    try:
        # check/remove header and empty lines
        if lines[0][0] == '-':
            lines = [line for line in lines[1:] if line]
        else:
            lines = [line for line in lines if line]
        data = [line.split() for line in lines]
        # remove non-amino acids
        data = [dat for dat in data if dat[2] == 'Amino']
        positions = [int(dat[0]) for dat in data]
        residues = [dat[3] for dat in data]
    except IndexError:
        error = 'E: unable to parse positions in file {}\n'.format(filename)
        sys.stderr.write(error)
        return False
    return [positions, residues]


def read_map_raw(indata):
    '''
    Read a map of tab-seperated alignment and data
    [0] for each line should be title, [1:] should be the data
    '''
    heads = list()
    seqs = list()
    data = list()
    for line in indata.split('\n'):
        if not line:
            continue
        cols = line.split('\t')
        if cols[0][:4] == 'ali:':
            heads.append(cols[0])
            seqs.append(''.join(cols[1:]))
        else:
            data.append(line)

    return (heads, seqs, data)


if __name__ == '__main__':
    if len(sys.argv) != 5:
        USAGE = ('Usage: {} <fasta alignment or data map> '.format(sys.argv[0]) +
                 '<protein name> <structure name> <ICM residue list>\n')
        sys.stderr.write(USAGE)
        sys.exit(1)
    if main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]) is False:
        sys.exit(2)
