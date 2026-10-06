#!/bin/bash

#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=128
#SBATCH --time=288:00:00
#SBATCH --mem=200g
#SBATCH --job-name=083_eggnog
#SBATCH --output=logs/083_eggnog.out
#SBATCH --error=logs/083_eggnog.err

module load devel/miniforge
conda activate eggnog_v2.1.12

cd /pfs/10/work/ho_yogau97-PiFerm3/PiFerm3_metagenomics/ 

#environment variables
BPMiFoDB=00_data/08_prokka_MiFoDB
BPcFMD=00_data/08_prokka_cFMD
BPdual=00_data/08_prokka_dual

BOUTMiFoDB=00_data/08_eggnog_MiFoDB
BOUTcFMD=00_data/08_eggnog_cFMD
BOUTdual=00_data/08_eggnog_dual

BOUTMiFoDBreport=02_report/08_eggnog_MiFoDB
BOUTcFMDreport=02_report/08_eggnog_cFMD
BOUTdualreport=02_report/08_eggnog_dual

#create directory
mkdir -p $BOUTMiFoDB
mkdir -p $BOUTcFMD
mkdir -p $BOUTdual

mkdir -p $BOUTMiFoDBreport
mkdir -p $BOUTcFMDreport
mkdir -p $BOUTdualreport

# add to path
export PATH=~/.conda/envs/eggnog_v2.1.12:~/.conda/envs/eggnog_v2.1.12/bin:"$PATH"

# download eggnog databases
export EGGNOG_DATA_DIR=$TMPDIR
download_eggnog_data.py -y

# run emapper per sample of MiFoDB
for FILE in $BPMiFoDB/*.faa; do 
 # sampleid
 MAG=$(basename "$FILE")
 MAG="${MAG/.faa/}" 
 
 if [ -f "$BOUTMiFoDB/$MAG.emapper.hits" ]; then #check if output files already generated
	 echo "already processed $MAG"
 else
 	 echo "start emapper for $MAG"

 	 emapper.py -i $FILE --itype proteins -m diamond --no_file_comments --cpu 128 --dbmem --tax_scope "prokaryota_broad" --temp_dir $TMPDIR --output_dir $BOUTMiFoDB -o $MAG # as in embl ebi pipeline
 	 echo "processed $MAG"
 fi
done

# run emapper per sample of cFMD 
for FILE in $BPcFMD/*.faa; do 
 # sampleid
 MAG=$(basename "$FILE")
 MAG="${MAG/.faa/}" 
 
 if [ -f "$BOUTcFMD/$MAG.emapper.hits" ]; then #check if output files already generated
	 echo "already processed $MAG"
 else
 	 echo "start emapper for $MAG"

 	 emapper.py -i $FILE --itype proteins -m diamond --no_file_comments --cpu 128 --dbmem --tax_scope "prokaryota_broad" --temp_dir $TMPDIR --output_dir $BOUTcFMD -o $MAG # as in embl ebi pipeline
 	 echo "processed $MAG"
 fi
done

# run emapper per sample of dual
for FILE in $BPdual/*.faa; do 
 # sampleid
 MAG=$(basename "$FILE")
 MAG="${MAG/.faa/}" 
 
 if [ -f "$BOUTdual/$MAG.emapper.hits" ]; then #check if output files already generated
	 echo "already processed $MAG"
 else
 	 echo "start emapper for $MAG"

 	 emapper.py -i $FILE --itype proteins -m diamond --no_file_comments --cpu 128 --dbmem --tax_scope "prokaryota_broad" --temp_dir $TMPDIR --output_dir $BOUTdual -o $MAG # as in embl ebi pipeline
 	 echo "processed $MAG"
 fi
done

# copy relevant files to report
cp $BOUTMiFoDB/*.annotations $BOUTMiFoDBreport/
cp $BOUTcFMD/*.annotations $BOUTcFMDreport/
cp $BOUTdual/*.annotations $BOUTdualreport/
