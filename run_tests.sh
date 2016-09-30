#!/bin/bash

PYTHONPATH=. py.test -vv tests/test_bioinfo.py --cov=bioinfo\
	                 tests/test_conservation.py --cov=conservation\
	                 tests/test_sequence_filter.py --cov=sequence_filter\
	                 tests/test_sort_by_blast.py --cov=sort_by_blast\
	                 tests/test_species.py --cov=species\
	                 tests/test_taxnode.py --cov=taxnode\
