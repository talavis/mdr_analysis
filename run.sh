#!/bin/bash

if [[ $# -ne 2 ]] ; then
    echo Usage: $0 sequencefile blastdbname
    exit 1
fi

# 1: input for BLAST
# 2: database location
BLASTIN=$1
DB=$2

PROJ=${BLASTIN%.*}
PROJ=${PROJ##*/}
QUERY=${BLASTIN##*/}
QUERYNAME=`head -1 ${BLASTIN} | cut -d '|' -f 2`
RUNDIR=run-${PROJ}-`date +%y%m%d-%H:%M`

if [[ ${0:0:1} == "/" ]] ; then
    SCRIPTDIR=${0%/*}
else
    SCRIPTDIR=`pwd`/${0%/*}
fi

mkdir ${RUNDIR}
cp ${BLASTIN} ${RUNDIR}
cd ${RUNDIR}

blastp -query ${QUERY} -db ${DB} -outfmt "6 sacc slen evalue pident length" -out ${PROJ}_blastp -num_threads `nproc` -num_alignments 1000
cut -f 1 ${PROJ}_blastp > ${PROJ}_accs
blastdbcmd -entry_batch ${PROJ}_accs -db ${DB} -out ${PROJ}_hits.fa
${SCRIPTDIR}/sequence_filter.py ${PROJ}_hits.fa ${QUERYNAME} > ${PROJ}_filtered
mafft-linsi --thread `nproc` ${PROJ}_filtered > ${PROJ}_filtered.mali
${SCRIPTDIR}/conservation.py ${PROJ}_filtered.mali ${QUERYNAME} > ${PROJ}_conservation

cd ..
