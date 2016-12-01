#!/usr/bin/env python3
'''
Map conserved residues to the active site.
Input:
1: Protein sequence
2: Structure code
3: ICM visualisation (gen_icm_visual)
4: ICM positions for active site
Output: Protein positions matching active site
'''

import sys

import bioinfo
import map_as_to_ali as maa


def main(prot_file, structcode, icmvis_file, icmpos_file):
    '''
    Map conserved residues to the active site.
    '''
    protseq = bioinfo.read_fasta(prot_file)
    if protseq is False:
        return False
    protseq = protseq[1][0]

    structseq = bioinfo.get_structseq(structcode)
    if structseq is False:
        return False
    structseq = structseq[1][0]

    conspos, consres = read_vis(icmvis_file)
    if conspos is False:
        return False

    icmpos = maa.read_icmdata(icmpos_file)
    if icmpos is False:
        return False
    icmpos = icmpos[0]

    data = [p for p in zip(conspos, consres) if p[0] in icmpos]
    # pos1 == 0
    pos = [d[0]-1 for d in data]
    res = [d[1] for d in data]

    prot_pos = maa.map_sequences(structseq, protseq, pos, res)
    if prot_pos is False:
        return False

    for i in range(len(prot_pos)):
        print('{}\t{}'.format(prot_pos[i],
                              res[i]))


def read_vis(visfile):
    '''
    Read the output from gen_icm_visual and
    return the positions and residues
    '''
    with open(visfile) as infile:
        raw = infile.read()
    line = [line for line in raw.split('\n') if line[:4] == 'cons'][0]
    rawpos = [pos.replace(',', '') for pos in line.split('^')[1:]]

    pos = [int(p[1:]) for p in rawpos]
    res = [p[0] for p in rawpos]
    if len(pos) == 0:
        err = 'E: No positions in ICM command file ({})\n'.format(visfile)
        sys.stderr.write(err)
        return False
    return pos, res


if __name__ == '__main__':
    if len(sys.argv) != 5:
        USAGE = ('Usage: {} '.format(sys.argv[0]) +
                 '<protein sequence> ' +
                 '<structure code> ' +
                 '<ICM visualisation> ' +
                 '<ICM active site positions>\n')
        sys.stderr.write(USAGE)
        sys.exit(1)

    if main(*sys.argv[1:]) is False:
        sys.exit(1)
