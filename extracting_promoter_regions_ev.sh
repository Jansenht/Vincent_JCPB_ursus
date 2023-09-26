# Download scripts from Github
git clone https://github.com/RimGubaev/extract_promoters.git
 
# NOTE: Script depends on bedtools and samtools to run. Might need to install these or load their modules if running on Kamiak. 
 
cd extract_promoters
 
idev # start interactive session on kamiak
module load bedtools
module load samtools
module load r
 
# Start script. It's interactive, so it'll prompt you for promoter length (I did 1000), genome fasta location, and gff file location
bash ./extract_promoters.sh
 
# I wrote a quick python script to rename the promoter fasta file to be a bit clearer and to include gene symbols.
# The output of this script is promoters.1kb.GCF_023065955.1_UrsArc1.0.fa
 