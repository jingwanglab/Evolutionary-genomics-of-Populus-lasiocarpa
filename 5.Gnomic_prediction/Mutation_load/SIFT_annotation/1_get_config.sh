#!?bin/bash
OutDir="/usr_storage/lzq/work/mutation/sift/sift4g"
SIFT4G="/usr_storage/lzq/software/sift4g/bin/sift4g"
echo -e "" > plas.sift.config.txt
echo -e "GENETIC_CODE_TABLE=1" >> plas.sift.config.txt
echo -e "GENETIC_CODE_TABLENAME=Standard" >> plas.sift.config.txt
echo -e "MITO_GENETIC_CODE_TABLE=11" >> plas.sift.config.txt
echo -e "MITO_GENETIC_CODE_TABLENAME=Plant Plastid Code" >> plas.sift.config.txt
echo -e "" >> plas.sift.config.txt
echo -e "PARENT_DIR=$OutDir" >> plas.sift.config.txt
echo -e "ORG=Populus_lasiocarpa" >> plas.sift.config.txt
echo -e "ORG_VERSION=v2" >> plas.sift.config.txt
echo -e "DBSNP_VCF_FILE=" >> plas.sift.config.txt
echo -e "" >> plas.sift.config.txt

echo -e "#Running SIFT 4G" >> plas.sift.config.txt
echo -e "SIFT4G_PATH=$SIFT4G" >> plas.sift.config.txt
echo -e "PROTEIN_DB=$OutDir/gene-annotation-src/uniprot_sprot.fasta" >> plas.sift.config.txt
echo -e "COMPUTER=GIS-KATNISS" >> plas.sift.config.txt
echo -e "" >> plas.sift.config.txt
echo -e "# Sub-directories, don't need to change" >> plas.sift.config.txt

echo -e "GENE_DOWNLOAD_DEST=gene-annotation-src" >> plas.sift.config.txt
echo -e "CHR_DOWNLOAD_DEST=chr-src" >> plas.sift.config.txt
echo -e "LOGFILE=Log.txt" >> plas.sift.config.txt
echo -e "ZLOGFILE=Log2.txt" >> plas.sift.config.txt
echo -e "FASTA_DIR=fasta" >> plas.sift.config.txt
echo -e "SUBST_DIR=subst" >> plas.sift.config.txt
echo -e "ALIGN_DIR=SIFT_alignments" >> plas.sift.config.txt
echo -e "SIFT_SCORE_DIR=SIFT_predictions" >> plas.sift.config.txt
echo -e "SINGLE_REC_BY_CHR_DIR=singleRecords" >> plas.sift.config.txt
echo -e "SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores" >> plas.sift.config.txt
echo -e "DBSNP_DIR=dbSNP" >> plas.sift.config.txt
echo -e "" >> plas.sift.config.txt
echo -e "# Doesn't need to change" >> plas.sift.config.txt
echo -e "FASTA_LOG=fasta.log" >> plas.sift.config.txt
echo -e "INVALID_LOG=invalid.log" >> plas.sift.config.txt
echo -e "PEPTIDE_LOG=peptide.log" >> plas.sift.config.txt
echo -e "ENS_PATTERN=ENS" >> plas.sift.config.txt
echo -e "SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord" >> plas.sift.config.txt

