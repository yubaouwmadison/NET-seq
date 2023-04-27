
#environment setup
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install R
conda install -c bioconda bowtie
#conda install -c bioconda cutadapt 
conda create -n cutadaptenv cutadapt
conda install -c bioconda samtools=1.9 --force-reinstall

#trimming (input=*.fastq.gz, output=*finaltrim.fastq)
conda activate cutadaptenv

cutadapt -q 20 -o  qualitytrim.fastq YBNET-WT-2_S34_L004_R1_001.fastq.gz  --cores=20   ##file name should be

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -o nodownstream.fastq  qualitytrim.fastq --cores=20

cutadapt -a CACAGTCAGTGTACATCGCTAAGTGACT -o noGAB.fastq   nodownstream.fastq --cores=20

cutadapt -a TCGCGTTGCGGGTCGACTCCGTGTACAT  -o nocontrol.fastq  noGAB.fastq --cores=20

cutadapt -g TCCGACGATCATTGATGGTGCCTACAG -o noadapter.fastq nocontrol.fastq --cores=20

cutadapt --max-n 0 -o  completetrim.fastq  noadapter.fastq --cores=20

cutadapt -m 10  -o  finaltrim.fastq completetrim.fastq --cores=20

rm qualitytrim.fastq nodownstream.fastq noGAB.fastq nocontrol.fastq noadapter.fastq completetrim.fastq 

conda deactivate

#bowtie alignment (input=*finaltrim.fastq, output=desiredregion.bed)
mkdir bowtie
cd bowtie

##copy MG1655v3.fasta to the folder /wt2/bowtie#

bowtie-build -f MG1655v3.fasta -q MG1655v3

cp finaltrim.fastq /mnt/scratch/yubao/netseq/wt2/bowtie

bowtie -x -q MG1655v3 finaltrim.fastq -S alignbowtie2.sam --best --strata -a -M 1 -v 2 --threads 20

samtools view -bSh -O BAM -o desiredregion2.bam  alignbowtie2.sam -@ 15

samtools view -bh -O BAM -o rrnE.bam -U non_rrnE.bam -L rrnE_Filter_RNAListMG1655v3.bed desiredregion2.bam -@ 15

bedtools bamtobed -i rrnE.bam > rrnE.bed


#make text file in R (input=desiredregion.bed, output=*reads.txt)

R #in the same folder start R

WT_BED<-read.table("rrnE.bed", sep="\t", header=FALSE)

WT_BEDTechnicalPlus<-WT_BED[which(WT_BED[,6]=="+"),]
WT_BEDTechnicalMinus<-WT_BED[which(WT_BED[,6]=="-"),]

WT_BEDTechnicalPlus[,3]<-as.integer(WT_BEDTechnicalPlus[,2]+15)
WT_BEDTechnicalMinus[,2]<-as.integer(WT_BEDTechnicalMinus[,3]-15)

WT_TrimBED<-rbind(WT_BEDTechnicalPlus, WT_BEDTechnicalMinus)

write.table(WT_TrimBED, file="rrnE_trim.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
rm(WT_TrimBED)

## copy MG1655v3BED_GenomeFile.genome to this folder

bedtools genomecov -i rrnE_trim.bed -d -strand - -g MG1655v3BED_GenomeFile.genome -5 > All3PrimePlus.txt
bedtools genomecov -i rrnE_trim.bed -d -strand + -g MG1655v3BED_GenomeFile.genome -5 > All3PrimeMinus.txt
bedtools genomecov -i rrnE_trim.bed -d -strand - -g MG1655v3BED_GenomeFile.genome > AllCoveragePlus.txt
bedtools genomecov -i rrnE_trim.bed -d -strand + -g MG1655v3BED_GenomeFile.genome > AllCoverageMinus.txt

mkdir rrn_txtfiles
mv All3PrimePlus.txt All3PrimeMinus.txt AllCoveragePlus.txt AllCoverageMinus.txt  rrn_txtfiles/
cd rrn_txtfiles/



#process with R 
R

AllCoveragePlus<-read.table("AllCoveragePlus.txt", header=FALSE, sep="\t")
AllCoverageMinus<-read.table("AllCoverageMinus.txt", header=FALSE, sep="\t")
All3PrimePlus<-read.table("All3PrimePlus.txt", header=FALSE, sep="\t")
All3PrimeMinus<-read.table("All3PrimeMinus.txt", header=FALSE, sep="\t")

Strand<-c(rep("+", 4641652), rep("-", 4641652)) 

Positions<-c(AllCoveragePlus[,2], AllCoverageMinus[,2])
Reads_Total<-c(AllCoveragePlus[,3], AllCoverageMinus[,3])  
Reads_3Prime<-c(All3PrimePlus[,3], All3PrimeMinus[,3])
Ratio_3Prime<-Reads_3Prime/Reads_Total
Ratio_3Prime[is.na(Ratio_3Prime)]<-0

DF_WT<-data.frame(Strand, Positions, Reads_Total, Reads_3Prime, Ratio_3Prime, stringsAsFactors=FALSE)    
write.table(DF_WT, file="readsandratio_rrnE.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#rrnE pause peak height
table<-read.table("readsandratio_rrnE.txt",header=TRUE, sep="\t")
Position_List<-c(4207532:4213234) #actual range:(4208146:4213159)
table<-table[which(table[,1]=="+"),]
plot(table[Position_List, c(2,4)], type="h", cex.axis=1.2, cex.lab=1.2, xlab="genomic position", ylab="normalized reads", ylim=c(0,10000))


