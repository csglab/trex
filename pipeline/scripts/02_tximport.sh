#!/usr/bin/bash

#SBATCH --job-name="txi"
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --err=logs/txi.%j.err
#SBATCH --output=outs/txi.%j.out
#SBATCH --mem-per-cpu=24G
#SBATCH --account=rrg-hsn
#SBATCH --array=0-30

# Input 

GROUPS=../example_inputs/groups.txt

# Modules

module load StdEnv/2020 r/4.2.1

# Build jobarray

GROUPS=()
i=0
while read p
do
  GROUPS[i]=$p
  i=$((i+1))
done < ${GROUPS}
GROUP=${GROUPS[${SLURM_ARRAY_TASK_ID}]}
OUTDIR=../results/${GROUP}/tximport

mkdir -p ${OUTDIR}

# Run tximport

echo "Processing sample ${GROUP}"
INDIR=../output/${GROUP}/salmon
MTDAT=../input/${GROUP}/metadata.csv

# Transcript counts
Rscript build-count-matrix.r ${INDIR} ${OUTDIR} ${MTDAT} ${GROUP}

# Gene counts
Rscript build-gene-count-matrix.r ${INDIR} ${OUTDIR} ${MTDAT} ${GROUP}