#import triplicate pause calling data, which containing redundant pause sites
wt2p<-read.table("/wt2_plus.txt",header=TRUE, sep="\t")
wt3p<-read.table("/wt3_plus.txt",header=TRUE, sep="\t")
wt4p<-read.table("/wt4_plus.txt",header=TRUE, sep="\t")

wt2m<-read.table("/wt2_minus.txt",header=TRUE, sep="\t")
wt3m<-read.table("/wt3_minus.txt",header=TRUE, sep="\t")
wt4m<-read.table("/wt4_minus.txt",header=TRUE, sep="\t")
#get unique intesect signals from triplicates

wt2p<-unique(wt2p)
wt2p<-wt2p[order(wt2p$Peak_List, decreasing=FALSE),]
wt3p<-unique(wt3p)
wt3p<-wt3p[order(wt3p$Peak_List, decreasing=FALSE),]
wt4p<-unique(wt4p)
wt4p<-wt4p[order(wt4p$Peak_List, decreasing=FALSE),]
wt2m<-unique(wt2m)
wt2m<-wt2m[order(wt2m$Peak_List, decreasing=FALSE),]
wt3m<-unique(wt3m)
wt3m<-wt3m[order(wt3m$Peak_List, decreasing=FALSE),]
wt4m<-unique(wt4m)
wt4m<-wt4m[order(wt4m$Peak_List, decreasing=FALSE),]

wtm<-unique(intersect(wt2m[,2], intersect(wt3m[,2], wt4m[,2])))
wtp<-unique(intersect(wt2p[,2], intersect(wt3p[,2], wt4p[,2])))
length(wtp)
length(wtm)

#import filter list containing tRNAs and sRNAs

filter_plus<-read.table("/filter_plus.txt",header=TRUE, sep="\t")
filter_minus<-read.table("/filter_minus.txt",header=TRUE, sep="\t")

#filter sites in filter list


plus<-wtp
minus<-wtm

site_plus<-c()
for (i in 1:length(plus)){
	for (j in 1:length(filter_plus[,1])){
		if(plus[i]>=filter_plus[j,1]&plus[i]<=filter_plus[j,2]){
			site_plus<-append(site_plus, plus[i])
			break
		}
	}
}

site_minus<-c()
for (i in 1:length(minus)){
	for (j in 1:length(filter_minus[,1])){
		if(minus[i]>=filter_minus[j,1]&minus[i]<=filter_minus[j,2]){
			site_minus<-append(site_minus, minus[i])
			break
		}
	}
}
length(site_plus)
length(site_minus)
wtp_f<-setdiff(wtp, site_plus)
wtm_f<-setdiff(wtm, site_minus)

#import list containing termination sites
tts_plus<-read.table("/tts_plus.txt", header=FALSE, sep="\t")
tts_minus<-read.table("/tts_minus.txt", header=FALSE, sep="\t")
tts_plus<-tts_plus[,1]
tts_minus<-tts_minus[,1]

#remove termination signals from pause signals
plus<-wtp_f
minus<-wtm_f

site_plus<-c()
for (i in 1:length(plus)){
	for (j in 1:length(tts_plus)){
		if(plus[i]==tts_plus[j]){
			site_plus<-append(site_plus, plus[i])
			break
		}
	}
}

site_minus<-c()
for (i in 1:length(minus)){
	for (j in 1:length(tts_minus)){
		if(minus[i]==tts_minus[j]){
			site_minus<-append(site_minus, minus[i])
			break
		}
	}
}
length(site_plus)
length(site_minus)
wtp_ff<-setdiff(wtp_f, site_plus)
wtm_ff<-setdiff(wtm_f, site_minus)
wtp_ff<-sort(wtp_ff, decreasing=FALSE)
wtm_ff<-sort(wtm_ff, decreasing=FALSE)

# remove redundants and sort by position for original pause data



strand<-rep("+",length(wtp_ff))
site<-wtp_ff

mean_3prime_read<-(wt2p[(which(wt2p[,2]%in%wtp_ff)),3]+wt3p[(which(wt3p[,2]%in%wtp_ff)),3]+wt4p[(which(wt4p[,2]%in%wtp_ff)),3])/3
mean_ratio<-(wt2p[(which(wt2p[,2]%in%wtp_ff)),4]+wt3p[(which(wt3p[,2]%in%wtp_ff)),4]+wt4p[(which(wt4p[,2]%in%wtp_ff)),4])/3

wtp_mean<-data.frame(strand, site, mean_3prime_read,mean_ratio)

strand<-rep("-",length(wtm_ff))
site<-wtm_ff
mean_3prime_read<-(wt2m[(which(wt2m[,2]%in%wtm_ff)),3]+wt3m[(which(wt3m[,2]%in%wtm_ff)),3]+wt4m[(which(wt4m[,2]%in%wtm_ff)),3])/3
mean_ratio<-(wt2m[(which(wt2m[,2]%in%wtm_ff)),4]+wt3m[(which(wt3m[,2]%in%wtm_ff)),4]+wt4m[(which(wt4m[,2]%in%wtm_ff)),4])/3
wtm_mean<-data.frame(strand, site, mean_3prime_read,mean_ratio)

write.table(wtp_mean, file="/wtp_mean.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(wtm_mean, file="/wtm_mean.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")


#add a 40 read confidence threshold
wtp_mean_40<-wtp_mean[which(wtp_mean[,3]>=40),]
wtm_mean_40<-wtm_mean[which(wtm_mean[,3]>=40),]
write.table(wtp_mean_40, file="/wtp_mean_40.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
write.table(wtm_mean_40, file="/wtm_mean_40.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")





