#!/bin/bash

TESTDIR=tests

PYTHONPATH=.:find_sequences py.test -vvv ${TESTDIR}/test_analyse_ali_groups.py --cov=analyse_ali_groups\
					 ${TESTDIR}/test_bioinfo.py --cov=bioinfo\
					 ${TESTDIR}/test_blast_by_name.py --cov=blast_by_name\
					 ${TESTDIR}/test_conservation.py --cov=conservation\
					 ${TESTDIR}/test_gen_icm_visual.py --cov=gen_icm_visual\
					 ${TESTDIR}/test_map_as_to_ali.py --cov=map_as_to_ali\
					 ${TESTDIR}/test_map_cons_to_ali.py --cov=map_cons_to_ali\
 					 ${TESTDIR}/test_map_icm_cons_as.py --cov=map_icm_cons_as\
					 ${TESTDIR}/test_select_conserved.py --cov=select_conserved\
					 ${TESTDIR}/test_sequence_filter.py --cov=sequence_filter\
					 ${TESTDIR}/test_seq_by_species.py --cov=seq_by_species\
					 ${TESTDIR}/test_seq_id.py --cov=seq_id\
					 ${TESTDIR}/test_sort_by_blast.py --cov=sort_by_blast\
					 ${TESTDIR}/test_species.py --cov=species\
					 ${TESTDIR}/test_taxnode.py --cov=taxnode\
					 ${TESTDIR}/find_sequences/test_select_accs_blastp.py --cov=select_accs_blastp

