#!/bin/bash

PYTHONPATH=..:../find_sequences py.test -vv test_bioinfo.py --cov=bioinfo\
	                                    test_conservation.py --cov=conservation\
	                                    test_sequence_filter.py --cov=sequence_filter\
	                                    test_sort_by_blast.py --cov=sort_by_blast\
	                                    test_species.py --cov=species\
	                                    test_taxnode.py --cov=taxnode\
                                            find_sequences/test_select_accs_blastp.py --cov=select_accs_blastp\
