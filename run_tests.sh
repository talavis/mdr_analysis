#!/bin/bash

TESTDIR=tests

PYTHONPATH=.:find_sequences py.test -vvv ${TESTDIR}/test_bioinfo.py --cov=bioinfo\
	                                 ${TESTDIR}/test_conservation.py --cov=conservation\
	                                 ${TESTDIR}/test_select_conserved.py --cov=select_conserved\
					 ${TESTDIR}/test_sequence_filter.py --cov=sequence_filter\
					 ${TESTDIR}/test_seq_id.py --cov=seq_id\
	                                 ${TESTDIR}/test_sort_by_blast.py --cov=sort_by_blast\
	                                 ${TESTDIR}/test_species.py --cov=species\
	                                 ${TESTDIR}/test_taxnode.py --cov=taxnode\
                                         ${TESTDIR}/find_sequences/test_select_accs_blastp.py --cov=select_accs_blastp\
