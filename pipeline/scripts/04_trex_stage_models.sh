#!/usr/bin/bash
#SBATCH --job-name="stage"
#SBATCH --err=stage.%j.err
#SBATCH --output=stage.%j.out
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --array=0-139

# Modules

module load StdEnv/2020 r/4.2.1

# Inputs

STAGE_PARS=params_stage.tsv
TREX_MAP=../references/gencode_GRCh38.p13/SUPPA2_splicemap/gencode.v37.primary_assembly.annotation_ASALL_strict.ioe # Splicing annotations including all events from SUPPA2

# Creating jobarray
params=()
i=0
while read p
do
  params[i]="$p"
  i=$((i+1))      
done < ${STAGE_PARS}

CANCER=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f1 )
EVENT=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f2 )
EVENT_COUNTS=../output/${CANCER}/trex/event_counts/${EVENT}.eventCounts.RData
OUT_DIR=../output/${CANCER}/trex/coefs_condition

mkdir -p ${OUT_DIR}

# Run TRex

Rscript trex_general-stage.r -s ${CANCER} \
                             -e ${EVENT} \
                             -c ${EVENT_COUNTS} \
                             -d ${OUT_DIR} \
                             -r ${TREX_MAP}