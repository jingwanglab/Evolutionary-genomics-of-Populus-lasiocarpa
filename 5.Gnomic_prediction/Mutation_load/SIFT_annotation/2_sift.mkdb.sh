sift="/usr_storage/lzq/software/sift4g/bin/sift4g"
gag="/usr_storage/lzq/software/GAG-master/gag.py"
fa="/usr_storage/lzq/work/mutation/snpeff/ZSP1912-L.LG.genome"
gff="/usr_storage/lzq/work/mutation/snpeff/genes.gff"
gffread="/usr_storage/lzq/software/gffread-0.12.7.Linux_x86_64/gffread"
python2 $gag --fasta $fa --gff $gff --fix_start_stop --out start_stop_fixed
cat start_stop_fixed/genome.comments.gff start_stop_fixed/genome.gff >start_stop_fixed/gene-annotation-src/plas.genes.gff3
cp start_stop_fixed/genome.proteins.fasta start_stop_fixed/gene-annotation-src/plas.pep.fasta
cp start_stop_fixed/genome.fasta start_stop_fixed/gene-annotation-src/plas.genome.fa
$gffread start_stop_fixed/gene-annotation-src/plas.genes.gff3 -T -o start_stop_fixed/gene-annotation-src/plas.genes.gtf
mkdb="/usr_storage/lzq/software/SIFT4G_Create_Genomic_DB-master/make-SIFT-db-all.pl"
cd /usr_storage/lzq/software/SIFT4G_Create_Genomic_DB-master/
####test the tool
###/usr/bin/perl make-SIFT-db-all.pl --config   test_files/homo_sapiens-test.txt 
perl make-SIFT-db-all.pl --config /usr_storage/lzq/work/mutation/sift/start_stop_fixed/plas.sift.config.txt
 fasta_formatter -i plas.genome.fa -w 60 -o plas.genome.new.fa
