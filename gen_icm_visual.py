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
    Generate the commands needed for visualisation of data on a structure in ICM
    pos 1 = 1
    '''
    protseq = bioinfo.read_fasta(protfile)[1][0]
    if protseq is False:
        return False
    structseq = bioinfo.get_structseq(structname)[1][0]
    if structseq is False:
        return False
    prot_pos = [int(num)-1 for num in open(posfile).read().split('\n') if len(num) > 0]
    prot_res = bioinfo.res_at_pos(protseq, prot_pos)

    struct_pos = maa.map_sequences(protseq, structseq, prot_pos, prot_res)
    if struct_pos is False:
        return False
    struct_res = bioinfo.res_at_pos(structseq, struct_pos)
    if struct_res is False:
        return False

    # make group of residues:
    icm_res = [make_icm_res(dat[0], dat[1], structname)
               for dat in zip(struct_res, struct_pos)]
    command = 'cons = {' + ' '.join(icm_res) + '}'
    print(command)

    
def make_icm_res(res, pos, structname):
    '''
    Return the formating for the residue in ICM
    pos 1 = 0
    '''
    return 'a_{}./^{}{}'.format(structname,
                                res,
                                pos+1)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        usage = 'Usage: {} <{}> <{}> <{}>\n'.format(sys.argv[0],
                                                    'Protein sequence file',
                                                    'PDB structure name',
                                                    'List of positions file')
        sys.exit(1)
    if main(sys.argv[1], sys.argv[2], sys.argv[3]) is False:
        sys.exit(1)
