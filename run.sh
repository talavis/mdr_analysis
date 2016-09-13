#!/bin/bash

SCRIPTDIR=.
${SCRIPTDIR}/sort_by_blast.py $1 $2 > ${1}_sorted
${SCRIPTDIR}/filter.py ${1}_sorted $3 > ${1}_filtered
mafft-linsi --thread `nproc` ${1}_filtered > ${1}_filtered.mali
