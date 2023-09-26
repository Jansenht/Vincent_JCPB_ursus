#!/bin/bash
#SBATCH --job-name=sam2bam
#SBATCH --partition=cas,kamiak
#SBATCH --output=/scratch/user/ellery.vincent/20220515_173202/1out.txt
#SBATCH --error=/scratch/user/ellery.vincent/20220515_173202/1err.txt
#SBATCH --workdir=/scratch/user/ellery.vincent/20220515_173202
#SBATCH --mem-per-cpu=100G
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job


# load samtools
module load samtools

# reference our path
path="/scratch/user/ellery.vincent/20220515_173202"


cd $path/

for i in *.sam
do
    #converting sam to bam
    samtools view -S -b $i > $path/$i.bam
    #Sort bam file by coordinate (takes .bam file and sorts it into a sorted bam file)
    samtools sort $path/$i.bam -o $path/$i.sorted.bam
    #compressing file to save space (optional)
    gzip $path/$i.bam
    gzip$path/$i.sorted.bam
done

echo "files converted to bam"