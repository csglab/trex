#!/usr/bin/bash
#SBATCH --job-name="counts"
#SBATCH --err=counts.%j.err
#SBATCH --output=counts.%j.out
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=1
#SBATCH --time=02:30:00
#SBATCH --array=0-90
#SBATCH --account=rrg-hsn

# Load modules

module load StdEnv/2020 r/4.2.1

# Define variables

COUNT_PARS=params_counts.mis
TREX_MAP=../references/gencode_GRCh38.p13/SUPPA2_splicemap/gencode.v37.primary_assembly.annotation_ASALL_strict.ioe # Splicing annotations including all events from SUPPA2
OUTFORMAT=tximeta # Output format of trex object

# Creating samples array
params=()
i=0
while read p
do
  params[i]=$p
  i=$((i+1))
done < ${COUNT_PARS}

CANCER=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f1 )
EVENT=$(echo ${params[${SLURM_ARRAY_TASK_ID}]} | cut -d "-" -f2 )

OUT_DIR=../output/${CANCER}/trex/event_counts
METADATA=../input/${CANCER}/metadata.csv
TRANSCRIPT_COUNTS=../output/${CANCER}/tximport/${CANCER}.txiobject.RData

mkdir -p ${OUT_DIR}

# Run TRex step 1

Rscript trex_general-event_counts.r -s ${CANCER} \
                                     -e ${EVENT} \
                                     -t ${OUTFORMAT} \
                                     -c ${TRANSCRIPT_COUNTS} \
                                     -m ${METADATA} \
                                     -r ${TREX_MAP} \
                                     -o ${OUT_DIR}