#!/bin/bash
#SBATCH --job-name=starmap
#SBATCH --partition=cas,kamiak
#SBATCH --output=/scratch/user/ellery.vincent/20220515_173202/1out.txt
#SBATCH --error=/scratch/user/ellery.vincent/20220515_173202/1err.txt
#SBATCH --workdir=/scratch/user/ellery.vincent/20220507_092042/star/trim_map_output
#SBATCH --mem-per-cpu=100G
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job



module load star


# grabs all files in the directory, displays how many files are in directory

FILES=(*.gz)
LEN=${#FILES[*]}

# loops through all files on a paired basis and runs STAR on paired reads

for (( i = 0; i < ${LEN}; i += 2 )); do
    echo "Running ${FILES[i]} and ${FILES[i+1]}"
    /opt/apps/star/2.7.6a/bin/STAR --genomeDir /data/lab/jansen/ursus_arctos_genome --runThreadN 20 --readFilesIn "${FILES[i]}" "${FILES[i+1]}" --outSAMstrandField intronMotif --readFilesCommand zcat --outMultimapperOrder Random --outSAMmultNmax 1 --outFilterMultimapNmax 10000 --outFileNamePrefix /scratch/user/ellery.vincent/20220515_173202/${FILES[i]}
    echo "* * * *"
done

echo "Script complete!"