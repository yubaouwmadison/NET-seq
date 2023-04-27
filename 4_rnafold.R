#import files with change analysis DOWN and UP
#import library Biostrings
library("Biostrings")
library("seqRFLP")
library("stringr")
MG1655_RawSeq<-readDNAStringSet("/MG1655v3.fasta", format="fasta")
MG1655Seq<-unlist(MG1655_RawSeq)
MG1655v3Raw<-toString(MG1655Seq)
MG1655v3DNA<-DNAString(MG1655v3Raw)
MG1655v3RNA<-RNAString(MG1655v3DNA)

table<-read.table("/table_ptest.txt",header=TRUE, sep="\t")
plus<-table[which(table[,1]=="+"&table[,10]=="DOWN"),2]
minus<-table[which(table[,1]=="-"&table[,10]=="DOWN"),2]

seq_list<-c()
for(i in 1:length(plus)){
 seq_tolist<-noquote(toString(MG1655v3RNA[(plus[i]-49):(plus[i]-10)]))
 seq_list<-append(seq_list, seq_tolist)
}
plus_seq<-seq_list

seq_list<-c()
for(i in 1:length(minus)){
 seq_tolist<-noquote(toString(reverse(complement(MG1655v3RNA[(minus[i]+10):(minus[i]+49)]))))
 seq_list<-append(seq_list, seq_tolist)
}
minus_seq<-seq_list

site<-plus
sequence<-plus_seq
df<-data.frame(site, sequence)
df$site<-sub("^","plus_",df$site)
temp<-dataframe2fas(df, file="/Users/yubao/Desktop/triple/full/4_rnafold/plus_40nt.fasta")

site<-minus
sequence<-minus_seq
df<-data.frame(site, sequence)
df$site<-sub("^","minus_",df$site)
temp<-dataframe2fas(df, file="/Users/yubao/Desktop/triple/full/4_rnafold/minus_40nt.fasta")

#enter terminal into the directory of plus_40nt.fasta
RNAfold < plus_40nt.fasta
RNAfold < minus_40nt.fasta

rm *.ps #remove the structure images
 command+s #save the result log of the terminal as *pred_plus.txt

#then repeat for the up changed pause signals line 16 to 47
plus<-table[which(table[,1]=="+"&table[,10]=="UP"),2]
minus<-table[which(table[,1]=="-"&table[,10]=="UP"),2]

# make tables from the previous results, for plus strand
rna<-read.table("Desktop/triple/full/4_rnafold/down_pred_plus.txt",header=FALSE, sep="\t")
x<-length(rna[,1])/3

site=c() #site of RNA at 3' end
rna_seq=c() #sequence of RNA in proximity of RNA exit channel -11 to -41
rnafold=c()
for (i in c(1:x)) {
	site_tolist<-rna[3*i-2,1]
	seq_tolist<-rna[3*i-1,1]
	fold_tolist<-rna[3*i,1]
	site<-append(site, site_tolist)
	rna_seq<-append(rna_seq, seq_tolist)
	rnafold<-append(rnafold, fold_tolist)
}

fold<-substr(rnafold, 1, 40)
energy<-as.numeric(substr(rnafold, 43,48))

site<-as.numeric(str_remove_all(site, "[>plus_]"))

#make better names for columns
plus_position<-site
exit_channel<-rna_seq
fold_prediction<-fold
free_energy<-energy

prediction<-data.frame(plus_position, exit_channel, fold_prediction, free_energy, stringsAsFactors=FALSE)
write.table(prediction, file="Desktop/triple/full/4_rnafold/down_fold_plus.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

# make tables from the previous results, for minus strand
rna<-read.table("Desktop/triple/full/4_rnafold/down_pred_minus.txt",header=FALSE, sep="\t")
x<-length(rna[,1])/3

site=c() #site of RNA at 3' end
rna_seq=c() #sequence of RNA in proximity of RNA exit channel -11 to -41
rnafold=c()
for (i in c(1:x)) {
	site_tolist<-rna[3*i-2,1]
	seq_tolist<-rna[3*i-1,1]
	fold_tolist<-rna[3*i,1]
	site<-append(site, site_tolist)
	rna_seq<-append(rna_seq, seq_tolist)
	rnafold<-append(rnafold, fold_tolist)
}

fold<-substr(rnafold, 1, 40)
energy<-as.numeric(substr(rnafold, 43,48))

site<-as.numeric(str_remove_all(site, "[>minus_]"))

#make better names for columns
plus_position<-site
exit_channel<-rna_seq
fold_prediction<-fold
free_energy<-energy

prediction<-data.frame(plus_position, exit_channel, fold_prediction, free_energy, stringsAsFactors=FALSE)
write.table(prediction, file="Desktop/triple/full/4_rnafold/down_fold_minus.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#according to RNAfold prediction, find pause hairpins: -11, -12, -13, -14, with at least 4 bp for stem loop and allowing at most 1 mismatch

plus_fold<-read.table("Desktop/triple/full/4_rnafold/down_fold_plus.txt", header=TRUE, sep="\t")
minus_fold<-read.table("Desktop/triple/full/4_rnafold/down_fold_minus.txt", header=TRUE, sep="\t")

p11<-plus_fold[which(substr(plus_fold[,3],37,40)=="))))"),]
p12<-plus_fold[which(substr(plus_fold[,3],36,40)=="))))."|substr(plus_fold[,3],35,40)=="))).)."|substr(plus_fold[,3],35,40)==")).))."|substr(plus_fold[,3],35,40)==").)))."),]
p13<-plus_fold[which(substr(plus_fold[,3],35,40)==")))).."|substr(plus_fold[,3],34,40)=="))).).."|substr(plus_fold[,3],34,40)==")).)).."|substr(plus_fold[,3],34,40)==").))).."),]
p14<-plus_fold[which(substr(plus_fold[,3],34,40)=="))))..."|substr(plus_fold[,3],33,40)=="))).)..."|substr(plus_fold[,3],33,40)==")).))..."|substr(plus_fold[,3],33,40)==").)))..."),]
plus_hairpin<-rbind(p11,p12,p13,p14)
length(p11[,1])
length(p12[,1])
length(p13[,1])
length(p14[,1])
plus_hairpin$hairpin_site<-c(rep("-11",164),rep("-12",228),rep("-13",166),rep("-14",122))
#numbers are respective PH numbers at certain positions

m11<-minus_fold[which(substr(minus_fold[,3],37,40)=="))))"),]
m12<-minus_fold[which(substr(minus_fold[,3],36,40)=="))))."|substr(minus_fold[,3],35,40)=="))).)."|substr(minus_fold[,3],35,40)==")).))."|substr(minus_fold[,3],35,40)==").)))."),]
m13<-minus_fold[which(substr(minus_fold[,3],35,40)==")))).."|substr(minus_fold[,3],34,40)=="))).).."|substr(minus_fold[,3],34,40)==")).)).."|substr(minus_fold[,3],34,40)==").))).."),]
m14<-minus_fold[which(substr(minus_fold[,3],34,40)=="))))..."|substr(minus_fold[,3],33,40)=="))).)..."|substr(minus_fold[,3],33,40)==")).))..."|substr(minus_fold[,3],33,40)==").)))..."),]
length(m11[,1])
length(m12[,1])
length(m13[,1])
length(m14[,1])
minus_hairpin<-rbind(m11,m12,m13,m14)

minus_hairpin$hairpin_site<-c(rep("-11",137),rep("-12",239),rep("-13",170),rep("-14",103))

write.table(plus_hairpin, file="Desktop/triple/full/4_rnafold/hpause_plus.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(minus_hairpin, file="Desktop/triple/full/4_rnafold/hpause_minus.txt", col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

