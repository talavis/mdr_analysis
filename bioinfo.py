#!/usr/bin/env python3
'''
Some commonly used bioinformatics funtions
'''

import re
import sys
import requests


def beautify_fasta(sequence, per_line=60):
    '''
    Limit the number of characters per line for a sequence
    '''
    i = per_line
    while i < len(sequence):
        sequence = sequence[:i] + '\n' + sequence[i:]
        i += per_line + 1

    return sequence


def get_species(header):
    '''
    Get the species of a sequence
    '''
    # NCBI
    if header[:2] == 'gi':
        return header[header.index('[')+1:header.index(']')]
    # UniProt
    if header[:2] in ['sp', 'tr']:
        return header[header.index('OS=')+3:
                      header.index('=', header.index('OS=')+3)-2].strip()

    sys.stderr.write('E: Unable to identify species({})'.format(header))

    return 'Unknown'


def get_structseq(pdb):
    '''
    Obtain the sequence of the protein structure from RCSB PDB
    '''
    url = ('http://www.rcsb.org/pdb/download/viewFastaFiles.do?' +
           'structureIdList={}&compressionType=uncompressed'.format(pdb))
    req = requests.get(url)
    raw = req.text
    if not raw or '<' in raw:
        error = 'E: could not retrieve sequence for structure {}\n'.format(pdb)
        sys.stderr.write(error)
        return False
    heads, seqs = read_fasta_raw(raw)

    return (heads, seqs)


def parse_icm_sel(selection):
    '''
    Parse an ICM selection
    '''
    if re.match(r'^a_.{4}\.[a-z]', selection) is None:
        return False

    raw = selection[selection.index('/')+1:]
    positions = [int(pos[2:]) for pos in raw.split(',')
                 if ':' not in pos]
    ranges = [pos for pos in raw.split(',')
              if ':' in pos]
    for i in enumerate(ranges):
        limits = ranges[i[0]].split(':')
        limits[0] = int(limits[0][2:])
        limits[1] = int(limits[1][2:])
        positions += list(range(limits[0], limits[1]+1))
    return sorted(positions)


def read_fasta(filename):
    '''
    Read a FASTA file.
    Return the headers and sequences as lists
    '''
    try:
        with open(filename) as infile:
            seqs = list()
            headers = list()
            for line in infile:
                if line[0] == '>':
                    headers.append(line[1:].strip())
                    seqs.append('')
                elif line.strip():
                    seqs[-1] += line.strip()
    except FileNotFoundError:
        error = 'E: Cannot find FASTA file: {}\n'.format(filename)
        sys.stderr.write(error)
        return False
    return (headers, seqs)


def read_fasta_raw(raw):
    '''
    Parse a FASTA formatted text string.
    Return the headers and sequences as lists
    '''
    seqs = list()
    headers = list()
    for line in raw.split('\n'):
        if not line:
            continue
        if line[0] == '>':
            headers.append(line[1:].strip())
            seqs.append('')
        elif line.strip():
            seqs[-1] += line.strip()

    return (headers, seqs)


def read_icm_res(icm_file):
    '''
    Parse ICM String(Res()) output
    '''
    sites = list()
    with open(icm_file) as infile:
        for line in infile.read().split('\n'):
            if not line:
                continue
            chain = re.compile(r'^a_.{4}\.[ab]')
            if '|' in line:
                parts = line.split('|')
            else:
                parts = [line]
            positions = list()
            for part in parts:
                if chain.match(part) is not None:
                    positions += parse_icm_sel(part)
            sites.append(sorted(positions))
    return sites


def res_at_pos(protseq, positions):
    '''
    Get a list of the residues at the positions in the protein
    pos 1 = 0
    '''
    if max(positions) > len(protseq) or min(positions) < 0:
        error = 'E: requested positions outside the protein sequence\n'
        sys.stderr.write(error)
        return False
    return [aa[1] for aa in enumerate(protseq) if aa[0] in positions]
