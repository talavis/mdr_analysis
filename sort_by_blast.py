#!/usr/bin/env python3

import sys
import bioinfo



def main() :
    if len(sys.argv) != 3 :
        sys.stderr.write('Usage: {} <BLAST tabular output> <FASTA file>\n'.format(sys.argv[0]))
        sys.exit(1)

    
        
if __name__ == '__main__' :
    main()
