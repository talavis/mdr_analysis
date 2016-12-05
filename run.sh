#!/bin/bash

if [[ $# -ne 4 ]] ; then
    echo Usage: $0 sequence.fa structure_code blastdb_path icm_pocket_file
    exit 1
fi

# 1: input for BLAST
# 2: database location
SEQUENCE=$1
STRUCTURE=$2
DB=$3
ICMPOCKET=$4

PROJ=${SEQUENCE%.*}
PROJ=${PROJ##*/}
QUERY=${SEQUENCE##*/}
POCKET=${ICMPOCKET##*/}
QUERYNAME=`head -1 ${SEQUENCE} | cut -d '|' -f 2`
RUNDIR=run-${PROJ}-`date +%y%m%d-%H:%M`

if [[ ${0:0:1} == "/" ]] ; then
    # absolute path
    SCRIPTDIR=${0%/*}
else
    # relative path
    SCRIPTDIR=`pwd`/${0%/*}
fi

mkdir ${RUNDIR}
cp ${SEQUENCE} ${RUNDIR}
cp ${ICMPOCKET} ${RUNDIR}
cd ${RUNDIR}

blastp -query ${QUERY} -db ${DB} -outfmt "6 sacc slen evalue pident length" -out ${PROJ}_blastp -num_threads `nproc` -num_alignments 1000 
cut -f 1 ${PROJ}_blastp > ${PROJ}_accs
blastdbcmd -entry_batch ${PROJ}_accs -db ${DB} -out ${PROJ}_hits.fa
${SCRIPTDIR}/blast_by_name.py ${PROJ}_hits.fa > ${PROJ}_genelim.fa
${SCRIPTDIR}/sequence_filter.py ${PROJ}_genelim.fa ${QUERYNAME} > ${PROJ}_filtered.fa
mafft-linsi --thread `nproc` ${PROJ}_filtered.fa > ${PROJ}_filtered.mali
${SCRIPTDIR}/conservation.py ${PROJ}_filtered.mali ${QUERYNAME} > ${PROJ}_conservation
${SCRIPTDIR}/select_conserved.py ${PROJ}_conservation 0.9 > ${PROJ}_cons90
cut -f 1 ${PROJ}_cons90 > ${PROJ}_cons90_pos
${SCRIPTDIR}/gen_icm_visual.py ${QUERY} ${STRUCTURE} ${PROJ}_cons90_pos > ${PROJ}_icm_cons_vis
${SCRIPTDIR}/map_icm_cons_as.py ${QUERY} ${STRUCTURE} ${PROJ}_icm_cons_vis ${POCKET} > ${PROJ}_active_site

cd ..
