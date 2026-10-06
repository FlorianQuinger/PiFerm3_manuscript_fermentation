#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --mem=120g
#SBATCH --job-name=082_prokka
#SBATCH --output=logs/082_prokka.out
#SBATCH --error=logs/082_prokka.err
cd /pfs/10/work/ho_yogau97-PiFerm3/PiFerm3_metagenomics

module load devel/miniforge
conda activate prokka_v1.14.6

#environment variables
DBdrepMiFoDB=00_data/08_db_drep_MiFoDB/dereplicated_genomes
DBdrepcFMD=00_data/08_db_drep_cFMD/dereplicated_genomes
DBdrepdual=00_data/08_db_drep_dual/dereplicated_genomes

BPMiFoDB=00_data/08_prokka_MiFoDB
BPcFMD=00_data/08_prokka_cFMD
BPdual=00_data/08_prokka_dual

#create directory
mkdir -p $BPMiFoDB
mkdir -p $BPcFMD
mkdir -p $BPdual

# loop over drepped files of MiFoDB and create faa

for FILE in $DBdrepMiFoDB/*; do
	MAG=$(basename "$FILE")
	MAG=${MAG/.fna/}
	prokka --outdir $BPMiFoDB --prefix $MAG --locustag $MAG --cpus 32 --force $FILE 
	echo "processed $MAG"
done

# loop over drepped files of cFMD and create faa

for FILE in $DBdrepcFMD/*; do
	MAG=$(basename "$FILE")
	MAG=${MAG/.fna/}
	prokka --outdir $BPcFMD --prefix $MAG --locustag $MAG --cpus 32 --force $FILE 
	echo "processed $MAG"
done

# loop over drepped files of dual DB and create faa

for FILE in $DBdrepdual/*; do
	MAG=$(basename "$FILE")
	MAG=${MAG/.fna/}
	prokka --outdir $BPdual --prefix $MAG --locustag $MAG --cpus 32 --force $FILE 
	echo "processed $MAG"
done