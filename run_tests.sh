#!/bin/bash

py.test -vv bioinfo.py conservation.py filter.py sort_by_blast.py species.py taxnode.py \
	--cov=bioinfo\
	--cov=conservation\
	--cov=filter\
	--cov=sort_by_blast\
	--cov=species\
	--cov=taxnode
