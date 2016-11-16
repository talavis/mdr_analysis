#!/usr/bin/env python3
'''
Generate the commands needed for visualisation of data on a structure in ICM

input: proteinsekvens, strukturnamn, positioner
output: kommandon för att visualisera i icm
'''

import sys

import bioinfo
import map_as_to_ali as maa


def main(protfile, structname, posfile):
    '''
    Generate the commands needed for visualisation of data
    on a structure in ICM
    pos 1 = 1
    '''
    protseq = bioinfo.read_fasta(protfile)
    if protseq is False or len(protseq) == 0:
        return False
    protseq = protseq[1][0]

    structseq = bioinfo.get_structseq(structname)
    if structseq is False:
        return False
    structseq = structseq[1][0]

    prot_pos = [int(num)-1 for num in
                open(posfile).read().split('\n')
                if len(num) > 0]
    prot_res = bioinfo.res_at_pos(protseq, prot_pos)

    struct_pos = maa.map_sequences(protseq, structseq, prot_pos, prot_res)
    if struct_pos is False:
        return False
    struct_res = bioinfo.res_at_pos(structseq, struct_pos)
    if struct_res is False:
        return False

    # read and convert the structure
    print('read pdb "{}"'.format(structname))
    print('convertObject a_{}. '.format(structname) +
          '1==1 no yes yes yes yes yes ' +
          '""+( 1==2 ? "water=tight ":"" )')

    # make group of residues:
    icm_res = [make_icm_res(dat[0], dat[1])
               for dat in zip(struct_res, struct_pos)]
    command = 'cons = a_{}./{}'.format(structname,
                                       ','.join(icm_res))
    print(command)

    # reset color to green
    print('color a_{}. green'.format(structname))

    # color all relevant residues
    print('color cons red')

    # activate graphical representation
    print('cool a_{}.'.format(structname))


def make_icm_res(res, pos):
    '''
    Return the formating for the residue in ICM
    pos 1 = 0
    '''
    return '^{}{}'.format(res, pos+1)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        USAGE = 'Usage: {} <{}> <{}> <{}>\n'.format(sys.argv[0],
                                                    'Protein sequence file',
                                                    'PDB structure name',
                                                    'List of positions file')
        sys.stderr.write(USAGE)
        sys.exit(1)
    if main(sys.argv[1], sys.argv[2], sys.argv[3]) is False:
        sys.exit(1)
