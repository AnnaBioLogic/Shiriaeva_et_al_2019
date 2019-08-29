#source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")

library(ShortRead)
setwd("D:/OneDrive/Results/Oligo transformation/HTS_new/1_2/PQ14 0N_longer105")



rm(list=ls())
new_spacers<-c(sread(readFasta("spacers.fasta")),sread(readFasta("maybe_spacers.fasta") ))
new_spacers<-DNAStringSet(new_spacers)



new_spacers<-gsub("R","N",new_spacers)
new_spacers<-gsub("Y","N",new_spacers)
new_spacers<-gsub("M","N",new_spacers)
new_spacers<-gsub("K","N",new_spacers)
new_spacers<-gsub("S","N",new_spacers)
new_spacers<-gsub("W","N",new_spacers)
new_spacers<-gsub("H","N",new_spacers)
new_spacers<-gsub("B","N",new_spacers)
new_spacers<-gsub("V","N",new_spacers)
new_spacers<-gsub("D","N",new_spacers)
new_spacers<-gsub("N","a",new_spacers)

length(new_spacers)#
################################

new_spacers<-DNAStringSet(new_spacers)
mychr1 <- sread(readFasta("pcdf_ec-cas12.fasta")) 
mychr_rev1<-reverseComplement(mychr1)


mychr2 <- sread(readFasta("bl21de3.fasta")) 
mychr_rev2<-reverseComplement(mychr2)


mychr3<-"NTCAAGCCCAATTTACTACTCGTTCTGGTGTTTCTCGTN" #TCAAGCCCAATTTACTACTCGTTCTGGTGTTTCTCGT is the 37nt oligonucleotide sequence
mychr3_rev<-reverseComplement(DNAString(mychr3))






new_spacers<-new_spacers[width(new_spacers)>=13]
length(new_spacers) #length of spacer that we will align
new_spacers<-DNAStringSet(new_spacers)

stat<-data.frame(table(width(new_spacers)))
read<-"smth"
Mychr1<-1
Mychr_rev1<-1
Mychr2<-1
Mychr_rev2<-1
Mychr3<-1
Mychr3_rev<-1




df<-data.frame(read, Mychr1,Mychr_rev1, Mychr2, Mychr_rev2,Mychr3, Mychr3_rev )
colnames(df)<-c("read", "mysearch1", "mysearch_rev1", "mysearch2", "mysearch_rev2", "mysearch3", "mysearch_rev3")
df<-df[-1,]

for (i in 1:nrow(stat)){
  w<-as.integer(as.character(stat[i,1]))
  reads<-new_spacers[width(new_spacers)==w]
  
  mypdict <- PDict(reads, max.mismatch=2)
  
  mysearch1 <- countPDict(mypdict, DNAString(as.character(mychr1)), max.mismatch=2) # Searches all dictionary entries against DNA1.
  mysearch_rev1<-countPDict(mypdict, DNAString(as.character(mychr_rev1)), max.mismatch=2) # Searches all dictionary entries against DNA1.
  
  
  mysearch2 <- countPDict(mypdict, DNAString(as.character(mychr2)), max.mismatch=2) # Searches all dictionary entries against DNA2.
  mysearch_rev2<-countPDict(mypdict, DNAString(as.character(mychr_rev2)), max.mismatch=2) # Searches all dictionary entries against DNA2.
  
  
  mysearch3 <- countPDict(mypdict, DNAString(as.character(mychr3)), max.mismatch=2) # Searches all dictionary entries against DNA3.
  mysearch_rev3<-countPDict(mypdict, DNAString(as.character(mychr3_rev)), max.mismatch=2) # Searches all dictionary entries against DNA3.
  
  
  df1<-data.frame(mysearch1)
  df2<-data.frame(mysearch_rev1)
  df3<-data.frame(mysearch2)
  df4<-data.frame(mysearch_rev2)
  df5<-data.frame(mysearch3)
  df6<-data.frame(mysearch_rev3)
  reads<-data.frame(as.character(reads))
  
  df_small<-cbind(reads, df1, df2, df3, df4, df5, df6)
  colnames(df_small)<-c("read", "mysearch1", "mysearch_rev1", "mysearch2", "mysearch_rev2", "mysearch3", "mysearch_rev3")
  df<-rbind(df, df_small)
  rm(df_small)
  rm(reads)
  rm(mysearch1)
  rm(mysearch_rev1)
  rm(mysearch_rev2)
  rm(mysearch_rev3)
  rm(mysearch2)
  rm(mysearch3)
  
}



sum<-df$mysearch1+df$mysearch_rev1+df$mysearch2 +df$mysearch_rev2+df$mysearch3+df$mysearch_rev3



nrow(df[sum==0,]) #not aligned spacers

sum_pl<-df$mysearch1+df$mysearch_rev1
sum_gen<-df$mysearch2+df$mysearch_rev2
sum_oligo<-df$mysearch3+df$mysearch_rev3



length(new_spacers[sum_pl>0&sum_gen>0&sum_oligo>0]) #aligned to all 3 templates
length(new_spacers[sum_pl>0&sum_gen>0&sum_oligo==0]) #aligned to plasmid and genome
length(new_spacers[sum_pl>0&sum_gen==0&sum_oligo>0]) #aligned to plasmid and oligo
length(new_spacers[sum_pl==0&sum_gen>0&sum_oligo>0]) #aligned to genome and oligo



length(new_spacers[sum_pl>0&sum_gen==0&sum_oligo==0]) #aligned to plasmid only
length(new_spacers[sum_pl==1&sum_gen==0&sum_oligo==0]) #unique on plasmid 


length(new_spacers[sum_pl==0&sum_gen>0&sum_oligo==0]) #aligned to genome only
length(new_spacers[sum_pl==0&sum_gen==1&sum_oligo==0]) #unique on genome 

length(new_spacers[sum_pl==0&sum_gen==0&sum_oligo>0]) #aligned to oligo only
length(new_spacers[sum_pl==0&sum_gen==0&sum_oligo==1]) #unique on oligo 








sp_ol_dir<-df$read[sum==1&df$mysearch3==1]
length(sp_ol_dir) #all oligo variants on dir strand


sp_ol_rev<-df$read[sum==1&df$mysearch_rev3==1]
length(sp_ol_rev) #all oligo variants on rev strand




#### align on dir oligo strand
reads<-DNAStringSet(sp_ol_dir)
mychr<-"NTCAAGCCCAATTTACTACTCGTTCTGGTGTTTCTCGTN"

F_amount <- countPDict(reads, DNAString(as.character(mychr)), max.mismatch=2)
table(F_amount)

reads<-reads[F_amount==1,]
mysearch <- matchPDict(reads, DNAString(as.character(mychr)), max.mismatch=2) # Searches all dictionary entries against chromosome.
Forward<-data.frame(mysearch)
Forward$count<-1

For<-aggregate(Forward$count, by=list(Forward$start, Forward$end, Forward$width), FUN=sum)
colnames(For)<-c("start", "end", "width", "count")

seq<-substr(rep(mychr, nrow(For)), For$start, For$end)
For<-cbind(For, seq) #table with spacers alignment to the direct strand of the oligo


###how much properly processed???
value<-For[(For$seq=="GCCCAATTTACTACTCGTTCTGGTGTTTCTCGT")&((For$end-For$start)==32),4]
value



For$orient<-"dir"





#### align on rev oligo strand
reads<-DNAStringSet(sp_ol_rev)
mychr<-("NACGAGAAACACCAGAACGAGTAGTAAATTGGGCTTGAN")

F_amount <- countPDict(reads, DNAString(as.character(mychr)), max.mismatch=2, fixed=FALSE)
table(F_amount)



reads<-reads[F_amount==1,]
mysearch_rev <- matchPDict(reads, DNAString(as.character(mychr)), max.mismatch=2, fixed=FALSE) # Searches all dictionary entries against chromosome.
Reverse<-data.frame(mysearch_rev)
Reverse$count<-1

Rev<-aggregate(Reverse$count, by=list(Reverse$start, Reverse$end, Reverse$width), FUN=sum)
colnames(Rev)<-c("start", "end", "width", "count")

seq<-substr(rep(mychr, nrow(Rev)), Rev$start, Rev$end)
Rev<-cbind(Rev, seq)#table with spacers alignment to the direct strand of the oligo


###how much properly processed???
value<-Rev[(Rev$seq=="ACGAGAAACACCAGAACGAGTAGTAAATTGGGC")&((Rev$end-Rev$start)==32),4]
value
Rev$orient<-"rev"





data<-rbind(For, Rev)
write.table(data, "Oligo_derived_spacers")