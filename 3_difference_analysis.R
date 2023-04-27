wt2<-read.table("/wt2.txt",header=TRUE, sep="\t")
wt3<-read.table("/wt3.txt",header=TRUE, sep="\t")
wt4<-read.table("/wt4.txt",header=TRUE, sep="\t")
at2<-read.table("/at2.txt",header=TRUE, sep="\t")
at3<-read.table("/at3.txt",header=TRUE, sep="\t")
at4<-read.table("/at4.txt",header=TRUE, sep="\t")

#normalizing with factors

wt2[,3]<-wt2[,3]*1.0562
wt2[,4]<-wt2[,4]*1.0562
wt3[,3]<-wt3[,3]*0.8560
wt3[,4]<-wt3[,4]*0.8560
wt4[,3]<-wt4[,3]*0.6977
wt4[,4]<-wt4[,4]*0.6977
at2[,3]<-at2[,3]*0.7406
at2[,4]<-at2[,4]*0.7406
at3[,3]<-at3[,3]*0.5191
at3[,4]<-at3[,4]*0.5191
at4[,3]<-at4[,3]*0.6488
at4[,4]<-at4[,4]*0.6488

wt2<-wt2[,-6]
wt3<-wt3[,-6]
wt4<-wt4[,-6]
at2<-at2[,-6]
at3<-at3[,-6]
at4<-at4[,-6]

write.table(wt2, file="/wt2.txt", col.names= TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(wt3, file="/wt3.txt", col.names= TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(wt4, file="/wt4.txt", col.names= TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(at2, file="/at2.txt", col.names= TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(at3, file="/at3.txt", col.names= TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(at4, file="/at4.txt", col.names= TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#first import wt2/wt3/wt4 at2/at3/at4 read values
#then import pause sites as plus and minus


cp<-
cbind(wt2[which(wt2[,1]=="+"&wt2[,2]%in%plus[,2]),c(3,4)],
wt3[which(wt3[,1]=="+"&wt3[,2]%in%plus[,2]),c(3,4)],
wt4[which(wt4[,1]=="+"&wt4[,2]%in%plus[,2]),c(3,4)],
at2[which(at2[,1]=="+"&at2[,2]%in%plus[,2]),c(3,4)],
at3[which(at3[,1]=="+"&at3[,2]%in%plus[,2]),c(3,4)],
at4[which(at4[,1]=="+"&at4[,2]%in%plus[,2]),c(3,4)]
) # reads total and reads 3prime from replicates, plus strand

cm<-
cbind(wt2[which(wt2[,1]=="-"&wt2[,2]%in%minus[,2]),c(3,4)],
wt3[which(wt3[,1]=="-"&wt3[,2]%in%minus[,2]),c(3,4)],
wt4[which(wt4[,1]=="-"&wt4[,2]%in%minus[,2]),c(3,4)],
at2[which(at2[,1]=="-"&at2[,2]%in%minus[,2]),c(3,4)],
at3[which(at3[,1]=="-"&at3[,2]%in%minus[,2]),c(3,4)],
at4[which(at4[,1]=="-"&at4[,2]%in%minus[,2]),c(3,4)]
) # reads total and reads 3prime from replicates, minus strand


#combine reads tatal and reads 3prime from triplicates
r1<-c()
r2<-c()
r3<-c()
r4<-c()
for (i in 1:5418)
{
		r1<-append(r1, cp[i,1]+cp[i,3]+cp[i,5])
		r2<-append(r2, cp[i,2]+cp[i,4]+cp[i,6])
		r3<-append(r3, cp[i,7]+cp[i,9]+cp[i,11])
		r4<-append(r4, cp[i,8]+cp[i,10]+cp[i,12])	
		
}
ratio_plus<-data.frame(plus[,c(1,2)],r1, r2, r3, r4)
ratio_plus<-ratio_plus[-which(ratio_plus[,6]==0|ratio_plus[,5]==0),] #remove zero read sites

r1<-c()
r2<-c()
r3<-c()
r4<-c()
for (i in 1:5286)
{
		r1<-append(r1, cm[i,1]+cm[i,3]+cm[i,5])
		r2<-append(r2, cm[i,2]+cm[i,4]+cm[i,6])
		r3<-append(r3, cm[i,7]+cm[i,9]+cm[i,11])
		r4<-append(r4, cm[i,8]+cm[i,10]+cm[i,12])	
		
}
ratio_minus<-data.frame(minus[,c(1,2)],r1, r2, r3, r4)
ratio_minus<-ratio_minus[-which(ratio_minus[,6]==0|ratio_minus[,5]==0),] #remove zero read sites

ratio<-rbind(ratio_plus, ratio_minus)

p_list<-c()
fold_list<-c()

for (i in 1: length(ratio[,1])){
	x1<-ratio[i,3]
	x2<-ratio[i,4]
	x3<-ratio[i,5]
	x4<-ratio[i,6]
	
	fold<-(x4/x3)/(x2/x1)
	fold_list<-append(fold_list, fold)
	
	p_read<-prop.test(x=c(x2,x4), n=c(x1,x3))$p.value
	p_list<-append(p_list, p_read)

}
q_list<-p.adjust(p_list, method = "BH", n = length(p_list))
logp<--log10(p_list)
logq<--log10(q_list)
log2fc<-log2(unlist(fold_list))
table_ptest<-data.frame(ratio[,c(1:6)], log2fc, logp, logq)  #proportion test for ratios

table_ptest$diff<-"NO"
table_ptest$diff[table_ptest$log2fc< -0.3040 & table_ptest$logq>4] <- "DOWN"
table_ptest$diff[table_ptest$log2fc> 0.4854 & table_ptest$logq>4] <- "UP"

#for plot aesthetics, hide some dots with extreme values
table_ptest[which(table_ptest[,8]=="Inf"),9]<-200
#some p values are too small to be displayed on plot, change p value to 200

table_limit<-table_ptest[-which(table_ptest[,9]>250|table_ptest[,9]<0.1),]

library(ggplot2)
library(MASS) 
library(scales)

p <- ggplot(data=table_limit, aes(x=log2fc, y=logq, col=diff)) + geom_point(alpha=0.2) + xlim(-6,4)  +xlab("Log2(fold change)") + ylab("-Log10(FDR)")+ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))

q<-p+geom_hline(yintercept=4, linetype="dashed", color = "black") +geom_vline(xintercept=c(-0.322, 0.485), linetype="dashed", color = "black")

p0<-subset(table_limit, site%in%c(258,2090066,2737665, 3950402))
p1<-subset(table_limit, site%in%c(83641,1799207,3852912))
p2<-rbind(p0,p1)
#r<-q + geom_point(data=p2, colour="black") + geom_text(data=p2, color="black",label=c("thrL","hisL","pheL","ilvL","leuL","pheM","ivbL"), hjust=1.2, size=2.5)
r<-q + geom_point(data=p2, colour="black")

mycolors<-c("cornflowerblue","coral","grey")
names(mycolors)<-c("DOWN","UP","NO")
s<-r+scale_color_manual(values=mycolors)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))











