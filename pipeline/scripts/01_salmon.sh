#!/usr/bin/bash
#SBATCH --job-name="salmon"
#SBATCH --cpus-per-task=10
#SBATCH --time=2:00:00
#SBATCH --err=salmon.%j.err
#SBATCH --output=salmon.%j.out
#SBATCH --mem=40G
#SBATCH --array=0-610

# User inputs 

REF=../references/GRCh38/GRCh38.d1.vd1.fa # Genome fasta
REFBED=../references/GRCh38/hg38_RefSeq.bed # Genome RefSeq annotations BED file 
FAIDX=../references/gencode_GRCh38.p13/salmon_genomeIndex # Directory with salmon index
SAMPLEINFO=../example_inputs/sample_info.txt # File with two colums: group and prefix of paired fastq files

# Configure environment

module --force purge 
module load nixpkgs/16.09
module load gcc/7.3.0
module load openmpi/3.1.4
module load salmon/1.3.0 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export global_start_time=$SECONDS

# Build jobarray

pargrid=()
i=0
while read u
do
	pargrid[i]=$u
	i=$((i+1))
done < ${SAMPLEINFO} 
pars=${pargrid[${SLURM_ARRAY_TASK_ID}]}

# Retrieve inputs

GROUP=$(echo ${pars} | cut -d " " -f1)
FILEID=$(echo ${pars} | cut -d " " -f2)

export TMPDIR=~/scratch/tmp/${GROUP}
mkdir -p ${TMPDIR}

### Quantify using Salmon
###################################################

echo "-----------------------------"
echo "> Quantify with salmon"

# Inputs

FQDIR=../output/${GROUP}/fastqs
FQ1=${FQDIR}/${FILEID}_R1.fastq # Fastq file with R1
FQ2=${FQDIR}/${FILEID}_R2.fastq # Fastq file with R2
OUTDIR=../output/${GROUP}/salmon
SALMONOUT=${OUTDIR}/${FILEID} # Directory with salmon outputs

mkdir -p ${OUTDIR}

# Run salmon

start_time=$SECONDS
salmon quant -p ${SLURM_CPUS_PER_TASK} \
             -l A \
             -1 ${FQ1} \
             -2 ${FQ2} \
             -o ${SALMONOUT} \
             -i ${FAIDX} \
             --validateMappings \
             --gcBias
             
elapsed=$(( SECONDS - start_time ))
echo "---> Finished salmon  in : $(date -ud "@$elapsed" +'%H hr %M min %S sec')"

echo "Inspect salmon output"
FILENAME=${SALMONOUT}/quant.sf
if [ -f "${FILENAME}" ];then
    if [ -s "${FILENAME}" ];then
        echo "Corresponding quant.sf file exists and is not empty. Removing intermediate files"
        rm ${FQ1} ${FQ2} ${BAMFILE}
        rm ${TMPDIR}/${FILEID}*
    else
        echo "Corresponding quant.sf file exists but is empty"
        echo "Reanalyze ${FILEID}"
        exit
    fi
else
    echo "File ${FILENAME} not exists"
    echo "Reanalyze ${FILEID}"
    exit
fi


global_elapsed=$(( SECONDS - global_start_time ))
echo "################################################################### "
echo "Completed in: $(date -ud "@$global_elapsed" +'%H hr %M min %S sec')"
echo "################################################################### "
