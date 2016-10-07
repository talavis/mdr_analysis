#!/usr/bin/python3

'''
Select accessions from a tabbed blastp result file
The fields are expected to be: "sseqid slen evalue pident length"
'''

import sys


def read_results(filename):
    '''
    Read a blastp file
    Returns a list of hits
    '''

    data = list()

    with open(filename) as infile:
        for line in infile:
            data.append(line.rstrip().split('\t'))

    return data


def filter_e(results, evalue=0.0):
    '''
    Filter the results based on E values
    Returns a set of accessions (column 0)
    '''
    matches = set()
    for res in results:
        if float(res[2]) <= evalue:
            matches.add(res[0])
    return matches


def main(filename, evalue=0.0):
    '''
    Read a blastp file and filter the results based on E values
    Prints the resulting set of accessions
    '''
    data = read_results(filename)
    accs = filter_e(data, evalue)

    for acc in accs:
        print(acc)


if __name__ == '__main__':
    if len(sys.argv) not in (2, 3):
        sys.stderr.write('Usage: {} '.format(sys.argv[0]) +
                         '<blastp file> ' +
                         '[e value limit]\n')
        sys.exit(1)
    if len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        main(sys.argv[1], float(sys.argv[2]))
