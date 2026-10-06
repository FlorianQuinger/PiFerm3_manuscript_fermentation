#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=64g
#SBATCH --job-name=084_db_creation
#SBATCH --output=logs/084_db_creation.out
#SBATCH --error=logs/084_db_creation.err
cd /pfs/10/work/ho_yogau97-PiFerm3/PiFerm3_metagenomics

module load devel/miniforge

#environment variables
DB=02_report/08_db

DBMiFoDB=$DB/ferment_MiFoDB
DBcFMD=$DB/ferment_cFMD
DBdual=$DB/ferment_dual

BPMiFoDB=00_data/08_prokka_MiFoDB
BPcFMD=00_data/08_prokka_cFMD
BPdual=00_data/08_prokka_dual

BOUTMiFoDB=00_data/08_eggnog_MiFoDB
BOUTcFMD=00_data/08_eggnog_cFMD
BOUTdual=00_data/08_eggnog_dual

#create directory
mkdir -p $DB
mkdir -p $DBMiFoDB
mkdir -p $DBcFMD
mkdir -p $DBdual

# copy all bins for MiFoDB

mkdir -p $DBMiFoDB/original_db # for faa files
mkdir -p $DBMiFoDB/eggNOG # for eggnog files

cp $BPMiFoDB/*.faa $DBMiFoDB/original_db
echo "copied all files from $BPMiFoDB"

for i in $BOUTMiFoDB/*.annotations; do 
	MAG=$(basename $i)
	MAG=${MAG/.emapper.annotations/}
	cp $i $DBMiFoDB/eggNOG/${MAG}_eggNOG.tsv
	echo "copied $i"
done

echo "finished creation of MiFoDB database"

# copy all bins for cFMD

mkdir -p $DBcFMD/original_db # for faa files
mkdir -p $DBcFMD/eggNOG # for eggnog files

cp $BPcFMD/*.faa $DBcFMD/original_db
echo "copied all files from $BPcFMD"

for i in $BOUTcFMD/*.annotations; do 
	MAG=$(basename $i)
	MAG=${MAG/.emapper.annotations/}
	cp $i $DBcFMD/eggNOG/${MAG}_eggNOG.tsv
	echo "copied $i"
done

echo "finished creation of cFMD database"

# copy all bins for dual

mkdir -p $DBdual/original_db # for faa files
mkdir -p $DBdual/eggNOG # for eggnog files

cp $BPdual/*.faa $DBdual/original_db
echo "copied all files from $BPdual"

for i in $BOUTdual/*.annotations; do 
	MAG=$(basename $i)
	MAG=${MAG/.emapper.annotations/}
	cp $i $DBdual/eggNOG/${MAG}_eggNOG.tsv
	echo "copied $i"
done

echo "finished creation of dual database"