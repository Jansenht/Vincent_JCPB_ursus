#!/bin/bash
#SBATCH --job-name=featurecounts3
#SBATCH --partition=cas,kamiak
#SBATCH --output=/scratch/user/ellery.vincent/20220630_142533/map/featurecounts/1out.txt
#SBATCH --error=/scratch/user/ellery.vincent/20220630_142533/map/featurecounts/1err.txt
#SBATCH --workdir=/scratch/user/ellery.vincent/20220630_142533/map/
#SBATCH --mem-per-cpu=70G
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --cpus-per-task=16

# load conda environment and activate feature counts within that environment

module load anaconda3
source activate featureCounts_env

echo "modules loaded succesfully"

# run featurecounts 

featureCounts *.sortedByCoord.out.bam -p -T 16 -g gene_id -a /data/lab/jansen/updated_genome/GCF_023065955.1_UrsArc1.0_genomic.gtf -o featureCountsgenelevel_output.txt 

echo "featurecounts file is generated"