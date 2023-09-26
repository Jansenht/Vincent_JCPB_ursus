#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --partition=kamiak
#SBATCH --output=/home/ellery.vincent/fastqc/out/1o
#SBATCH --error=/home/ellery.vincent/fastqc/out/1e
#SBATCH --workdir=/data/lab/jansen/RNA_raw_data
#SBATCH --mem-per-cpu=75G
#SBATCH --nodes=1                ### Number of nodes
#SBATCH --ntasks=1               ### Number of tasks per array job
#SBATCH --time=01-10:00:00

module load fastqc/0.11.9

cd /data/lab/jansen/RNA_raw_data

fileList=($(find `pwd` -type f -name '*.fastq.gz')) 

for i in ${fileList[@]}
    do

    fastqc --outdir /home/ellery.vincent/fastqc/out --noextract --nogroup $i
    
    done

    