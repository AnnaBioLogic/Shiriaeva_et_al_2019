source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("Biostrings")

library(ShortRead)
library(Biostrings)


rm(list=ls())

#Strain - KD675
#Biological replicate 1 - for this replicate 2x75 paired-read sequencing

###This is how this library is supposed to appear (N - N11 and N9 unique molecular identifiers (UMI) conjugated to 5' and 3' adapters; seq - insert): 
#AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATCTCAGNNNNNNNNNNNseqNNNNNNNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACgtggccATCTCGTATGCCGTCTTCTGCTTG
#index primer RPI20 (index - GGCCAC) (See Supplementary table 24)
#The forward read sequence is supposed to start with 4-letter barcode TCAG (part of i115 5' adapter, see supplementary table 23) and be followed by an N11 UMI, an insert, an N9 UMI, and 3' adapter sequence (starting with TGGAATTCTCGGGTG....)
Barcode_4nt<-as.character("TCAG")

folder_sample<-as.character("D:/I_F fragments/Final fragments analysis/May 2018/Replica 1")#specify pathway to a folder with data for your sample
setwd(folder_sample) 



gc()

############################   2*75

setwd(folder_sample)
setwd("2_75")
F_75<-as.character("KD675_FragSeq_Rep1_HTS1_2x75_Forward_reads") #insert the name of the forward reads
R_75<-as.character("KD675_FragSeq_Rep1_HTS1_2x75_Reverse_reads") #insert the name of the reverse reads


###forward reads
reads_F=readFastq(F_75)
length(reads_F) #how many reads
width_F<-as.numeric(width(reads_F))
stat_F<-as.data.frame(table(width_F))
rm(reads_F)


##reverse reads
gc()
reads_R=readFastq(R_75)
width_R<-as.numeric(width(reads_R))
stat_R<-as.data.frame(table(width_R))
rm(reads_R)


##we will analyze read quality by groups of reads with the same read length. Therefore to keep forward and reverse reads in the same order we will trimm each pair of forward and reverse reads so that each pair has the same read length
####min width for each pair
min_width<-pmin(width_F, width_R)
rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt", "min_width", "F_75", "R_75")))
gc()

##forward reads
reads_F=readFastq(F_75)
reads_F<-narrow( reads_F, 1, min_width)

stat<-as.data.frame(table(min_width))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$min_width<-as.numeric(as.character(stat$min_width))

PQcutoff <- 20 #we will substitute bases with PhredQuality less than this to N

For<-c()
for (i in 1:nrow(stat)){
  width<-stat[i,1]
  reads<-reads_F[width(reads_F)==width]
  
  l<-length(reads)
  
  if (l<=100000) {
    reads_F1<-reads
    seq_R1 <- sread(reads_F1) 
    q <- PhredQuality(quality(quality(reads_F1))) 
    MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
    
    positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
    N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
    N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet")
    rm(reads_F1)
    seq_R <- replaceLetterAt(seq_R1, positions, N_tr)
    For<-append(For, seq_R)  
  }
  else {
    whole<-l%/%100000
    rest<-l%%100000
    for (b in 1:(whole+1)) {
      if (b<(whole+1)) {
        reads_F1<-reads[(b*100000-99999):(b*100000)]
        seq_R1 <- sread(reads_F1) 
        q <- PhredQuality(quality(quality(reads_F1))) 
        MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
        
        positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
        N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
        N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet")
        rm(reads_F1)
        seq_R <- replaceLetterAt(seq_R1, positions, N_tr)
        For<-append(For, seq_R)      }
      
      else {
        reads_F1<-reads[(b*100000-99999):length(reads)]
        seq_R1 <- sread(reads_F1) 
        q <- PhredQuality(quality(quality(reads_F1))) 
        MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
        
        positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
        N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
        N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet")
        rm(reads_F1)
        seq_R <- replaceLetterAt(seq_R1, positions, N_tr)
        For<-append(For, seq_R)
      }
      
    }
    
  }
  
}
rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt", "min_width", "For", "R_75")))
writeFasta(For, "Forward_trimmed.fasta")  


#reverse reads
gc()
rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt", "min_width", "R_75")))


reads_R<-readFastq(R_75)
reads_R<-narrow( reads_R, 1, min_width)
stat<-as.data.frame(table(min_width))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$min_width<-as.numeric(as.character(stat$min_width))



PQcutoff <- 20 #we will substitute bases with PhredQuality less than this to N
Rev<-c()
for (i in 1:nrow(stat)){
  width<-stat[i,1]
  reads<-reads_R[width(reads_R)==width]
  
  l<-length(reads)
  
  if (l<=100000) {
    reads_F1<-reads
    seq_R1 <- sread(reads_F1) 
    q <- PhredQuality(quality(quality(reads_F1))) 
    MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
    
    positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
    N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
    N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet")
    rm(reads_F1)
    seq_R <- replaceLetterAt(seq_R1, positions, N_tr) 
    Rev<-append(Rev, seq_R)  
  }
  else {
    whole<-l%/%100000
    rest<-l%%100000
    
    for (b in 1:(whole+1)) {
      if (b<(whole+1)) {
        reads_F1<-reads[(b*100000-99999):(b*100000)]
        seq_R1 <- sread(reads_F1) 
        q <- PhredQuality(quality(quality(reads_F1))) 
        MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
        
        positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
        N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
        N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet")
        rm(reads_F1)
        seq_R <- replaceLetterAt(seq_R1, positions, N_tr) 
        Rev<-append(Rev, seq_R)}
      
      
      else {
        reads_F1<-reads[(b*100000-99999):length(reads)]
        seq_R1 <- sread(reads_F1) 
        q <- PhredQuality(quality(quality(reads_F1))) 
        MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
        
        positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
        N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
        N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet")
        rm(reads_F1)
        seq_R <- replaceLetterAt(seq_R1, positions, N_tr) 
        Rev<-append(Rev, seq_R)
      }
      
    }
    
  }
  
}


writeFasta(Rev, "Reverse_trimmed.fasta")  





rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt")))
gc()
For<-sread(readFasta("Forward_trimmed.fasta"))
Rev<-sread(readFasta("Reverse_trimmed.fasta"))
Rev<-reverseComplement(Rev)


#finding overlaps between forward and reverse reads
myalign <- pairwiseAlignment( Rev, For, type = "local", substitutionMatrix = NULL, fuzzyMatrix = NULL, 
                              gapOpening = -10, gapExtension = -4, scoreOnly = FALSE)

myalign1<-myalign[width(pattern(myalign))>=20]
length(myalign1)
For_20<-For[width(pattern(myalign))>=20]
Rev_20<-Rev[width(pattern(myalign))>=20]

####Combine overlapping reads to get the entire read sequence (starting with 4 nt barcode and ending with N9 UMI sequence)
Combined_insert<-paste(
  substr(For_20,
         1,end(subject(myalign[width(pattern(myalign))>=20]))),
  substr(Rev_20,           
         end(pattern(myalign[width(pattern(myalign))>=20]))+1,
         width(Rev_20)),
  sep=""
)


plot(table(width(Combined_insert)))



log<-width(Combined_insert)>=40
table(log)



Combined_insert<-Combined_insert[log]
writeFasta(DNAStringSet(Combined_insert), "Combined_2_75.fasta")


rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt")))
setwd(folder_sample)
setwd("2_75")
reads<-readFasta("Combined_2_75.fasta")
reads<-DNAStringSet(sread(reads))




###which contain 4nt barcode?

first_four<-substr(reads, 1, 4)
hits=vcountPattern(Barcode_4nt,first_four,
                   max.mismatch=1,with.indels=FALSE)
log<-hits==1
table(log)


reads<-reads[log]
reads<-substr(reads, 5, width(reads)) ##remove 5' barcode



###reads 16-100
log<-(width(reads)>=36)&(width(reads)<=120) ##for consistency between all samples we will analyze fragments 16-100 nt
table(log)
reads<-reads[log]



###we have a lot of reads for this sample so we will keep only the reads with very high quality (all positions with PhredQuality>20,no Ns)
###wchich contain N? keep 0Ns
hits=vcountPattern("N",reads,
                   max.mismatch=0,with.indels=FALSE)
log<-hits==0
table(log)
reads<-reads[log]

###Reads with identical sequence of the insert and molecular identifiers (N11 and N9) are overamplified. We will substitute number of these reads with 1.
reads<-table(reads)
reads<-data.frame(reads)
reads<-as.character(reads$reads)
length(reads)

plot(table(width(reads)))

reads1<-substr(reads, 12, width(reads)-9)


setwd(folder_sample) #make sure that your folder contains "kd675_rc.fasta" reference genome and sequences of plasmids "prsf_csy.fasta" and "pcas.fasta"
writeFasta(DNAStringSet(reads1), "reads_for_alignment.fasta")

plot(table(width(reads1)), main="Read length")




########################################################################
#########################################align reads to the genome and plasmids

rm(list=setdiff(ls(), c("folder_sample")))



insert<-readFasta("reads_for_alignment.fasta")
insert<-sread(insert)


######align to plasmids/genome
#subject for mapping

mychr1<-sread(readFasta("prsf_csy.fasta")) 
mychr1_rev<-reverseComplement(mychr1)

mychr2<-sread(readFasta("pcas.fasta")) 
mychr2_rev<-reverseComplement(mychr2)


mychr3 <- sread(readFasta("kd675_rc.fasta")) # Creates sample chromosome.
mychr3_rev<-reverseComplement(mychr3) #creates rev sample chromosome



insert<-gsub("R","N",insert)
insert<-gsub("Y","N",insert)
insert<-gsub("M","N",insert)
insert<-gsub("K","N",insert)
insert<-gsub("S","N",insert)
insert<-gsub("W","N",insert)
insert<-gsub("H","N",insert)
insert<-gsub("B","N",insert)
insert<-gsub("V","N",insert)
insert<-gsub("D","N",insert)
insert<-gsub("N", "a", insert)
insert<-DNAStringSet(insert)






width<-width(insert)
tab<-data.frame(table(width))



seq<-c("-")
F_amount1<-c("-")
  F_amount1_rev<-c("-")
F_amount2<-c("-")
  F_amount2_rev<-c("-")
  F_amount3<-  c("-")  
F_amount3_rev<-c("-")


table<-data.frame(seq, F_amount1,     F_amount1_rev, F_amount2,     F_amount2_rev, F_amount3,    
                  F_amount3_rev)
table<-table[-1,]

for (i in 1:nrow(tab)){
  w<-as.numeric(as.character((tab[i,1])))
  
  set<-insert[width(insert)==w]
  mypdict <- PDict(set, max.mismatch=2)
  
  F_amount1 <- countPDict(mypdict, DNAString(as.character(mychr1)), max.mismatch=2)
  F_amount1_rev <- countPDict(mypdict, DNAString(as.character(mychr1_rev)), max.mismatch=2)
  
  F_amount2 <- countPDict(mypdict, DNAString(as.character(mychr2)), max.mismatch=2)
  F_amount2_rev <- countPDict(mypdict, DNAString(as.character(mychr2_rev)), max.mismatch=2)
  
  F_amount3 <- countPDict(mypdict, DNAString(as.character(mychr3)), max.mismatch=2)
  F_amount3_rev <- countPDict(mypdict, DNAString(as.character(mychr3_rev)), max.mismatch=2)
  
  df<-data.frame(set, F_amount1, F_amount1_rev, F_amount2, F_amount2_rev, F_amount3, F_amount3_rev )
  colnames(df)<-colnames(table)
  table<-rbind(table, df)
  rm(df)
  rm(set)
  rm(mypdict)
  
  }




table$sum<-table$F_amount1+ table$F_amount1_rev+ table$F_amount2+ table$F_amount2_rev+ table$F_amount3+ table$F_amount3_rev #how many alignments for each read
table_0<-table[table$sum==0,]
nrow(table_0) #how many not aligned reads
reads_NA<-DNAStringSet(as.character(table_0$seq))
writeFasta(reads_NA, "not_aligned_reads.fasta")




table_genome<-table[(table$F_amount3>0|table$F_amount3_rev>0),]
table_genome_unique<-table_genome[table_genome$sum==1,]
nrow(table_genome) #how many aligned to the genome
nrow(table_genome_unique) #how many uniquely aligned to the genome

reads_dir_genome<-as.character(table_genome_unique$seq[table_genome_unique$F_amount3==1])
reads_rev_genome<-as.character(table_genome_unique$seq[table_genome_unique$F_amount3_rev==1])


writeFasta(DNAStringSet(reads_dir_genome), "Unique_reads_genome_direct.fasta")
writeFasta(DNAStringSet(reads_rev_genome), "Unique_reads_genome_reverse.fasta")




table_plasmid1<-table[(table$F_amount1>0|table$F_amount1_rev>0),]
table_plasmid1_unique<-table_plasmid1[table_plasmid1$sum==1,]
nrow(table_plasmid1) #how many aligned to the plasmid 1
nrow(table_plasmid1_unique) #how many uniquely aligned to the plasmid 1


reads_dir_plasmid1<-as.character(table_plasmid1_unique$seq[table_plasmid1_unique$F_amount1==1])
reads_rev_plasmid1<-as.character(table_plasmid1_unique$seq[table_plasmid1_unique$F_amount1_rev==1])


writeFasta(DNAStringSet(reads_dir_plasmid1), "Unique_reads_prsf_csy_direct.fasta")
writeFasta(DNAStringSet(reads_rev_plasmid1), "Unique_reads_prsf_csy_reverse.fasta")




table_plasmid2<-table[(table$F_amount2>0|table$F_amount2_rev>0),]
table_plasmid2_unique<-table_plasmid2[table_plasmid2$sum==1,]
nrow(table_plasmid2)#how many aligned to the plasmid 2
nrow(table_plasmid2_unique) #how many uniquely aligned to the plasmid 2

reads_dir_plasmid2<-as.character(table_plasmid2_unique$seq[table_plasmid2_unique$F_amount2==1])
reads_rev_plasmid2<-as.character(table_plasmid2_unique$seq[table_plasmid2_unique$F_amount2_rev==1])

writeFasta(DNAStringSet(reads_dir_plasmid2), "Unique_reads_pCas_direct.fasta")
writeFasta(DNAStringSet(reads_rev_plasmid2), "Unique_reads_pCas_reverse.fasta")

rm(list=ls())




##Determine position of alignment for each uniquely aligned read
####align unique genome reads on dir strand
mychr <- sread(readFasta("kd675_rc.fasta"))
mychr_rev<-reverseComplement(mychr)


insert_un<-sread(readFasta("Unique_reads_genome_direct.fasta"))


width<-width(insert_un)
tab<-data.frame(table(width))
MM<-2 #mismatches allowed
w<-as.numeric(as.character((tab[1,1])))
set<-insert_un[width(insert_un)==w]
mypdict <- PDict(set, max.mismatch=MM)

F_amount <- countPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM)

mysearch <- matchPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM) # Searches all dictionary entries against chromosome.

Forward<-data.frame(mysearch)


log<-F_amount==1
reads<-set[log]
Forward<-cbind(Forward, reads)

for (i in 2:nrow(tab)){
  w<-as.numeric(as.character((tab[i,1])))
  
  
  set<-insert_un[width(insert_un)==w]
  mypdict <- PDict(set, max.mismatch=MM)
  
  F_amount <- countPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM)
  mysearch <- matchPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
  
  Forward1<-data.frame(mysearch)

  log<-F_amount==1
  reads<-set[log]
  Forward1<-cbind(Forward1, reads)
  
 
  Forward<-rbind(Forward, Forward1)
}

#coordinates for extraction
first<-Forward$start
second<-Forward$end

log<-((first-20)>=1)&((second+20)<=width(mychr))
Forward<-Forward[log,]
first<-first[log]
second<-second[log]


M<-rep(mychr, nrow(Forward))
M<-DNAStringSet(M)
reads_dir<-narrow(M, first-20, second+20)


w<-width(reads_dir)

Forward$reads<-as.character(narrow(reads_dir, 21, w-20) )#the sequence of each insert
Forward$read_20_20<-as.character(reads_dir) #the sequence of each insert plus 20 bp upstream and downstream in the chromosome

Forward$strand<-"dir"
Forward<-Forward[,-c(1:2)] #this is a data frame with unique reads aligned to the dir strand of KD403, each read is shown separately so there could be several reads that are exactly the same in the data frame

#write.table(Forward, "Forward_alignment_table")



##rev
#coordinates for extraction

insert_un<-sread(readFasta("Unique_reads_genome_reverse.fasta"))


width<-width(insert_un)
tab<-data.frame(table(width))
MM<-2 #mismatches allowed
w<-as.numeric(as.character((tab[1,1])))
set<-insert_un[width(insert_un)==w]
mypdict <- PDict(set, max.mismatch=MM)


F_amount_rev <- countPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM)
mysearch_rev<-matchPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM) # Searches all dictionary entries against chromosome.


Reverse<-data.frame(mysearch_rev)



log<-F_amount_rev==1
reads<-set[log]
Reverse<-cbind(Reverse, reads)


for (i in 2:nrow(tab)){
  w<-as.numeric(as.character((tab[i,1])))
  
  
  set<-insert_un[width(insert_un)==w]
  mypdict <- PDict(set, max.mismatch=MM)
  
 
  F_amount_rev <- countPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM)
  
   mysearch_rev<-matchPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
  
  
  Reverse1<-data.frame(mysearch_rev)
  
  
  log<-F_amount_rev==1
  reads<-set[log]
  Reverse1<-cbind(Reverse1, reads)
  
  
  
  Reverse<-rbind(Reverse, Reverse1)
}


first<-Reverse$start
second<-Reverse$end

log<-((first-20)>=1)&((second+20)<=width(mychr))
Reverse<-Reverse[log,]
first<-first[log]
second<-second[log]


M<-rep(mychr_rev, nrow(Reverse))
M<-DNAStringSet(M)
reads_rev<-narrow(M, first-20, second+20)

w<-width(reads_rev)

Reverse$reads<-as.character(narrow(reads_rev, 21, w-20))#the sequence of each insert
Reverse$read_20_20<-as.character(reads_rev) #the sequence of each insert plus 20 bp upstream and downstream in the chromosome


#the current coordinates for reverse reads in the Reverse table are for mychr_rev sequence (so that the first position in mychr_rev is the last position in mychr)
#let's calculate the coordinates relative to mychr

first<-Reverse$start
second<-Reverse$end

whole_chr<-width(mychr)
first_real<-whole_chr-second+1
second_real<-whole_chr-first+1
first<-first_real
second<-second_real
Reverse$start<-first
Reverse$end<-second
Reverse$strand<-"rev"
Reverse<-Reverse[,-c(1:2)]  #this is a data frame with unique reads aligned to the rev strand of KD403, each read is shown separately so there could be several reads that are exactly the same in the data frame; the start coordinate is now given relative to the mychr sequence (not mychr_rev) so it is actually the coordinate of reverse read 3' end



#write.table(Reverse, "Reverse_alignment_table")

Both_strands<-rbind(Forward, Reverse)
colnames(Both_strands)<-c("start", "end", "width", "read_sequence", "(20nt_up_chromosome)+READ+(20nt_dw_chromosome)", "strand")
write.table(Both_strands, "Alignment_table_both_strands") #this is the table with all unique reads mapped to dir or rev strand


