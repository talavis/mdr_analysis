#!/bin/bash

py.test -vv bioinfo.py --cov=bioinfo\
	    conservation.py --cov=conservation\
	    filter.py --cov=filter\
	    sort_by_blast.py --cov=sort_by_blast\
	    species.py --cov=species\
	    taxnode.py --cov=taxnode\
