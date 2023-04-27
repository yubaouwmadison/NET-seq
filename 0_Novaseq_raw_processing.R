
#environment setup
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install R
conda install -c bioconda bowtie
conda install -c bioconda cutadapt
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

bowtie -x -q MG1655v3 finaltrim.fastq -S alignbowtie.sam --threads 20

samtools view -bh -O BAM -o qualityalignment.bam -U rejected.bam  -q 20 alignbowtie.sam

##copy Latest_FIlter_RNAListMG1655v3.bed to the folder /wt2/bowtie#

samtools view -bh -O BAM -o filterout.bam -U desiredregion.bam -L Latest_Filter_RNAListMG1655v3.bed qualityalignment.bam

bedtools bamtobed -i desiredregion.bam > desiredregion.bed

#make text file in R (input=desiredregion.bed, output=*reads.txt)

R #in the same folder start R

WT_BED<-read.table("desiredregion.bed", sep="\t", header=FALSE)

WT_BEDTechnicalPlus<-WT_BED[which(WT_BED[,6]=="+"),]
WT_BEDTechnicalMinus<-WT_BED[which(WT_BED[,6]=="-"),]

WT_BEDTechnicalPlus[,3]<-as.integer(WT_BEDTechnicalPlus[,2]+15)
WT_BEDTechnicalMinus[,2]<-as.integer(WT_BEDTechnicalMinus[,3]-15)

WT_TrimBED<-rbind(WT_BEDTechnicalPlus, WT_BEDTechnicalMinus)

write.table(WT_TrimBED, file="trim.bed", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

## copy MG1655v3BED_GenomeFile.genome to this folder

bedtools genomecov -i trim.bed -d -strand - -g MG1655v3BED_GenomeFile.genome -5 > All3PrimePlus.txt
bedtools genomecov -i trim.bed -d -strand + -g MG1655v3BED_GenomeFile.genome -5 > All3PrimeMinus.txt
bedtools genomecov -i trim.bed -d -strand - -g MG1655v3BED_GenomeFile.genome > AllCoveragePlus.txt
bedtools genomecov -i trim.bed -d -strand + -g MG1655v3BED_GenomeFile.genome > AllCoverageMinus.txt

mv All3PrimePlus.txt All3PrimeMinus.txt AllCoveragePlus.txt AllCoverageMinus.txt  txtfiles/

#process with R

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
write.table(DF_WT, file="readsandratio.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#count Bsubtilis spike-in reads
bbmap.sh in=finaltrim.fastq outm=junk.fastq outu=ecoli_unmapped.fastq ref=MG1655v3.fasta
bbmap.sh in=ecoli_unmapped.fastq sam=1.3 out=bsub_mapped.sam ambiguous=toss ref=Bsubtilisv3.fasta perfectmode=t

#scale 3 prime end reads with normalization factor

wt2<-read.table("wt2_ReadsandRatio.txt", header=TRUE, sep="\t")
wt2$norm_reads<-as.integer(0.8588*wt2$Reads_3Prime)    ##(pick the right factor for each dataset)
write.table(wt2, file="wt2_norm.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")


