#!/usr/bin/bash
#SBATCH --job-name="condition"
#SBATCH --err=logs/condition.%j.err
#SBATCH --output=outs/condition.%j.out
#SBATCH --account=rrg-hsn
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --time=3:00:00
#SBATCH --array=0-153

# Load modules

module load StdEnv/2020 r/4.2.1

# Define variables

CONDITION_PARS=params_condition.tsv
TREX_MAP=../references/gencode_GRCh38.p13/SUPPA2_splicemap/gencode.v37.primary_assembly.annotation_ASALL_strict.ioe # Splicing annotations including all events from SUPPA2

params=()
i=0
while read p
do
  params[i]="$p"
  i=$((i+1))      
done < ${CONDITION_PARS}

CANCER=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f1 )
EVENT=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f2 )
EVENT_COUNTS=../output/${CANCER}/trex/event_counts/${EVENT}.eventCounts.RData
OUT_DIR=../output/${CANCER}/trex/coefs_condition

mkdir -p ${OUT_DIR}

# Run TRex 

Rscript trex_general-condition.r -s ${CANCER} \
                                 -e ${EVENT} \
                                 -c ${EVENT_COUNTS} \
                                 -d ${OUT_DIR} \
                                 -r ${TREX_MAP}

