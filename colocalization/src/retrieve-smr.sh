#!/usr/bin/env bash

TISSUE=${1}
PVAL=${2}

SMR='/data/CARD/projects/singlecell_humanbrain/coloc/smr-1.3.1'
OUTDIR="${PWD}/data/eQTL"
TMPDIR="/lscratch/${SLURM_JOB_ID}"
eQTL_DIR='/data/CARD/projects/omicSynth/SMR_omics/eQTLs'
FILESTEM="2020-05-26-${TISSUE}"
BESD="${eQTL_DIR}/${FILESTEM}/${FILESTEM}-${CHR}-SMR-besd"

retrieve_smr() {
    local chr=${1}
    local pval=${2}
    local besd="${eQTL_DIR}/${FILESTEM}/${FILESTEM}-${chr}-SMR-besd"
    echo "retrieving SMR for ${TISSUE} chr ${chr}"

    # print ALL eQTL results regardless of p-value
    ${SMR} \
        --descriptive-cis \
        --beqtl-summary \
        ${besd} \
        --query ${pval} \
        --out ${TISSUE}_${chr}
    
    if [ ! -f "${OUTDIR}/${TISSUE}.tsv" ]; then
        cat ${TISSUE}_${chr}.txt > ${OUTDIR}/${TISSUE}.tsv
    else
        awk 'NR > 1' ${TISSUE}_${chr}.txt >> ${OUTDIR}/${TISSUE}.signif.tsv
    fi
        rm ${TISSUE}_${chr}.txt
}


# RUN
mkdir -p ${OUTDIR}
mkdir -p ${TMPDIR} && cd ${TMPDIR}

echo "Tissue type ${TISSUE}"

for CHR in $(seq 1 22); do
    retrieve_smr ${CHR} ${PVAL}
done

echo 'gzipping...'
gzip ${OUTDIR}/${TISSUE}.tsv

echo "all done"