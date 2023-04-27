#For minus strand

input<-c("~/1_raw_normalized_data/wt2.txt", "~/1_raw_normalized_data/wt3.txt", "~/1_raw_normalized_data/wt4.txt")

output<-c("~/1_raw_normalized_data/wt2_minus.txt", "~/1_raw_normalized_data/wt3_minus.txt", "~/1_raw_normalized_data/wt4_minus.txt")


Gene_Pause_Minus<-function(input, output){

for (m in 1:length(input)){

Gene_Table<-read.table("~/1_raw_normalized_data/tu_minus.txt", header=FALSE, sep="\t")
Read_Table<-read.table(input[m], head=TRUE, sep="\t") 

Minus_Table<-Read_Table[which(Read_Table[,1]=="-"),]

gene_list<-as.character(Gene_Table[,2])
gene_end<-Gene_Table[,1]
gene_start<-Gene_Table[,2]
FirstRound_3PrimeReads<-Minus_Table[,6]
Ratio_Raw<-Minus_Table[,5]




Peak_List<-c()
Reads3Prime_List<-c()
DistancetoStart_List<-c()
Quality_List<-c()
GeneLength_List<-c()
Gene_Name<-c()
Ratio_List<-c()

for (i in 1:4){
  if (i==1){
    SecondRound_3PrimeReads<-FirstRound_3PrimeReads
    for (i in 1:length(gene_list)){
      positions<-c(gene_end[i]:gene_start[i])
      for(k in 1:length(positions)){
        Reads_List<-FirstRound_3PrimeReads[(positions[k]-100):(positions[k]+99)]
        Reads_NoZeros<-Reads_List[which(Reads_List>0)]
        if (length(Reads_NoZeros)>20&&FirstRound_3PrimeReads[positions[k]]!=0){
          Critical_Value<-mean(Reads_NoZeros)+4*sd(Reads_NoZeros)
          if (FirstRound_3PrimeReads[positions[k]]>Critical_Value){
            Peak_List<-append(Peak_List, positions[k])
            Reads3Prime_List<-append(Reads3Prime_List, FirstRound_3PrimeReads[positions[k]])
            DistancetoStart_List<-append(DistancetoStart_List, (gene_start[i]-positions[k]))
            GeneLength_List<-append(GeneLength_List, length(positions))
	    Quality_List<-append(Quality_List, (FirstRound_3PrimeReads[positions[k]])/Critical_Value)
	    Gene_Name<-append(Gene_Name, gene_list[i])
	    Ratio_List<-append(Ratio_List, Ratio_Raw[positions[k]])
            SecondRound_3PrimeReads[positions[k]]<-0
          }
        }
      }
    }
  }
  else if (i==2){
    ThirdRound_3PrimeReads<-SecondRound_3PrimeReads
    for (i in 1:length(gene_list)){
      positions<-c(gene_end[i]:gene_start[i])
      for(k in 1:length(positions)){
        Reads_List<-SecondRound_3PrimeReads[(positions[k]-100):(positions[k]+99)]
        Reads_NoZeros<-Reads_List[which(Reads_List>0)]
        if (length(Reads_NoZeros)>20&&SecondRound_3PrimeReads[positions[k]]!=0){
          Critical_Value<-mean(Reads_NoZeros)+4*sd(Reads_NoZeros)
          if (SecondRound_3PrimeReads[positions[k]]>Critical_Value){
            Peak_List<-append(Peak_List, positions[k])
            Reads3Prime_List<-append(Reads3Prime_List, FirstRound_3PrimeReads[positions[k]])
            DistancetoStart_List<-append(DistancetoStart_List, (gene_start[i]-positions[k]))
            GeneLength_List<-append(GeneLength_List, length(positions))
	    Quality_List<-append(Quality_List, (FirstRound_3PrimeReads[positions[k]])/Critical_Value)
	    Gene_Name<-append(Gene_Name, gene_list[i])
	    Ratio_List<-append(Ratio_List, Ratio_Raw[positions[k]])
            ThirdRound_3PrimeReads[positions[k]]<-0
          }
        }
      }
    }
  }
  else if (i==3){
    FourthRound_3PrimeReads<-ThirdRound_3PrimeReads
    for (i in 1:length(gene_list)){
      positions<-c(gene_end[i]:gene_start[i])
      for(k in 1:length(positions)){
        Reads_List<-ThirdRound_3PrimeReads[(positions[k]-100):(positions[k]+99)]
        Reads_NoZeros<-Reads_List[which(Reads_List>0)]
        if (length(Reads_NoZeros)>20&&ThirdRound_3PrimeReads[positions[k]]!=0){
          Critical_Value<-mean(Reads_NoZeros)+4*sd(Reads_NoZeros)
          if (ThirdRound_3PrimeReads[positions[k]]>Critical_Value){
            Peak_List<-append(Peak_List, positions[k])
            Reads3Prime_List<-append(Reads3Prime_List, FirstRound_3PrimeReads[positions[k]])
            DistancetoStart_List<-append(DistancetoStart_List, (gene_start[i]-positions[k]))
            GeneLength_List<-append(GeneLength_List, length(positions))
	    Quality_List<-append(Quality_List, (FirstRound_3PrimeReads[positions[k]])/Critical_Value)
	    Gene_Name<-append(Gene_Name, gene_list[i])
	    Ratio_List<-append(Ratio_List, Ratio_Raw[positions[k]])
            FourthRound_3PrimeReads[positions[k]]<-0
          }
        }
      }
    }
  }
  else if (i==4){
    for (i in 1:length(gene_list)){
      positions<-c(gene_end[i]:gene_start[i])
      for(k in 1:length(positions)){
        Reads_List<-FourthRound_3PrimeReads[(positions[k]-100):(positions[k]+99)]
        Reads_NoZeros<-Reads_List[which(Reads_List>0)]
        if (length(Reads_NoZeros)>20&&FourthRound_3PrimeReads[positions[k]]!=0){
          Critical_Value<-mean(Reads_NoZeros)+4*sd(Reads_NoZeros)
          if (FourthRound_3PrimeReads[positions[k]]>Critical_Value){
            Peak_List<-append(Peak_List, positions[k])
            Reads3Prime_List<-append(Reads3Prime_List, FirstRound_3PrimeReads[positions[k]])
            GeneLength_List<-append(GeneLength_List, length(positions))
            DistancetoStart_List<-append(DistancetoStart_List, (gene_start[i]-positions[k]))
	    Quality_List<-append(Quality_List, (FirstRound_3PrimeReads[positions[k]])/Critical_Value)
	    Ratio_List<-append(Ratio_List, Ratio_Raw[positions[k]])
	    Gene_Name<-append(Gene_Name, gene_list[i])
          }
        }
      }
    }
  }
}


Strand<-rep("-", length(Peak_List))

#Minus_DF_Window<-data.frame(Strand, Peak_List, DistancetoStart_List, Reads3Prime_List, Ratio_List, Quality_List, GeneLength_List, Gene_Name, stringsAsFactors=FALSE)
Minus_DF_Window<-data.frame(Strand, Peak_List, Reads3Prime_List, Ratio_List, stringsAsFactors=FALSE)
write.table(Minus_DF_Window, file=output[m], col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

}
}


Gene_Pause_Minus(input, output)

