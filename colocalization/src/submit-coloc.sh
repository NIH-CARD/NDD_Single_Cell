#!/usr/bin/env bash

module load R/4


for diagnosis in PD_Nalls_no23 PD_Nalls AD_Bellenguez LBD_Chia ALS_vanRheenen; do
    for tissue in Cerebellum-EUR Cortex-AFR Cortex-EUR Hippocampus-EUR Spinalcord-EUR Basalganglia-EUR; do
        echo $tissue $diagnosis
        Rscript ./run_coloc.R $tissue $diagnosis
    done
done

