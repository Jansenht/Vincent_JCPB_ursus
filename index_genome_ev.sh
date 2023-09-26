#!/bin/bash
#SBATCH --partition=cas,kamiak         ### Partition
#SBATCH --job-name=index2              ### Job name
#SBATCH --output=/data/lab/jansen/updated_genome/index/1out.txt      ### File in which to store job output
#SBATCH --error=/data/lab/jansen/updated_genome/index/1err.txt       ### File in which to store job error messages
#SBATCH --workdir=/data/lab/jansen/updated_genome
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node (MPI processes)
#SBATCH --mem=100G
#SBATCH --cpus-per-task=16
module load star

echo "Module loaded"

path="data/lab/jansen/updated_genome"

STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /data/lab/jansen/updated_genome/index --genomeFastaFiles /data/lab/jansen/updated_genome/GCF_023065955.1_UrsArc1.0_genomic.fna --sjdbGTFfile /data/lab/jansen/updated_genome/GCF_023065955.1_UrsArc1.0_genomic.gff

echo "Index genome generated"
