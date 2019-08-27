source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")


library(ShortRead)


rm(list=ls())

#Strain - KD403
#Biological replicate 3 - for this replicate 1x150 single-read sequencing was performed and 2x150 paired-read sequencing

###This is how this library is supposed to appear (N - N11 and N9 unique molecular identifiers (UMI) conjugated to 5' and 3' adapters; seq - insert): 
#AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATCAGTCNNNNNNNNNNNseqNNNNNNNNNTGGAATTCTCGGGTGCCAAGGAACTCCAGTCACggctacATCTCGTATGCCGTCTTCTGCTTG
#index primer RPI11 (index - GTAGCC) (See Supplementary table 24)
#The forward read sequence is supposed to start with 4-letter barcode AGTC (part of i114 5' adapter, see supplementary table 23) and be followed by an N11 UMI, an insert, an N9 UMI, and 3' adapter sequence (starting with TGGAATTCTCGGGTG....)
Barcode_4nt<-as.character("AGTC")



####1*150
folder_sample<-as.character("D:/FragSeq/I-E/3.3")#specify pathway to a folder with data for your sample
setwd(folder_sample)  
setwd("1_150") #specify pathway to a folder with 1x150 fastq file
fastq_1_150<-as.character("KD403_FragSeq_Rep3_HTS1_1x150")




reads_F=readFastq(fastq_1_150)
a<-length(reads_F)
a #number of reads 


#dividing big files into smaller portions (on my computer I can't analyze more than 15000000 reads)
reads_F<-reads_F[15000001:a]  

########################################################
#########################################   1.  Quality control

width_F<-as.numeric(width(reads_F))
stat<-as.data.frame(table(width_F))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$width_F<-as.numeric(as.character(stat$width_F))

PQcutoff <- 20 #we will substitute bases with PhredQuality less than this to N
rm(width_F)



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
        For<-append(For, seq_R)        }
      
      
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

rm(list=setdiff(ls(), c("For", "folder_sample", "fastq_1_150", "Barcode_4nt")))
seq<-For
rm(list=setdiff(ls(), c("seq", "folder_sample", "fastq_1_150", "Barcode_4nt")))



########################################################
#########################################   2.   Trimming 3' adapters
###wchich contain TGGAATTCTCGGGTG? 
hits_3prime=vcountPattern("TGGAATTCTCGGGTG",seq,
                          max.mismatch=1,with.indels=FALSE)
log<-hits_3prime==1
table(log)

####looking for RPI position, 1MM allowed
seq<-seq[log]
l<-length(seq)
div<-l/200000
div
div<-ceiling(div)
div

start<-c()
for(i in 1:(div-2)){
  min<-i*200000-199999
  max<-i*200000
  fastq<-seq[min:max]
  start_RPI<-startIndex(vmatchPattern("TGGAATTCTCGGGTG",fastq, max.mismatch=1,with.indels=FALSE))
  start_RPI<-as.vector(unlist(start_RPI))
  start<-append(start, start_RPI)
}
min<-(div-1)*200000-199999
max<-l
fastq<-seq[min:max]
start_RPI<-startIndex(vmatchPattern("TGGAATTCTCGGGTG",fastq, max.mismatch=1,with.indels=FALSE))
start_RPI<-as.vector(unlist(start_RPI))
start<-append(start, start_RPI) #this is the final vector with the start coordinate of RPI adapter (the border btw N9 and RPI)
rm("fastq")
length(start)

reads_trimmed<-narrow(x=seq, start = 1, end = (start-1)) 
plot(table(width(reads_trimmed)))



########################################################
########################################## 3.  Getting rid of adapter dimers and reads with too short inserts
reads_F<-reads_trimmed[width(reads_trimmed)>=40]  #keep the reads with an indert>=16
plot(table(width(reads_F)))
length(reads_F)
writeFasta(reads_F, "Forward_trimmed1.fasta")  #write sequences with an insert into a Fasta file (each read starts with 4nt barcode followed by N11, insert and N9)




###############################################################################
#####repeat steps 1-3 for the rest of reads 

rm(list=setdiff(ls(), c("folder_sample", "fastq_1_150", "Barcode_4nt")))
gc()
reads_F=readFastq(fastq_1_150)
reads_F<-reads_F[1:1500000]

########################################################
#########################################   1.  Quality control
width_F<-as.numeric(width(reads_F))
stat<-as.data.frame(table(width_F))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$width_F<-as.numeric(as.character(stat$width_F))

PQcutoff <- 20 #we will substitute bases with PhredQuality less than this to N
rm(width_F)

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

rm(list=setdiff(ls(), c("For","folder_sample","Barcode_4nt")))
seq<-For
rm(list=setdiff(ls(), c("seq", "folder_sample", "Barcode_4nt")))



########################################################
#########################################   2.   Trimming 3' adapters
###wchich contain TGGAATTCTCGGGTG?
hits_3prime=vcountPattern("TGGAATTCTCGGGTG",seq,
                          max.mismatch=1,with.indels=FALSE)
log<-hits_3prime==1
table(log)

####looking for RPI position, 1MM
seq<-seq[log]
l<-length(seq)
div<-l/200000
div
div<-ceiling(div)
div

start<-c()
for(i in 1:(div-2)){
  min<-i*200000-199999
  max<-i*200000
  fastq<-seq[min:max]
  start_RPI<-startIndex(vmatchPattern("TGGAATTCTCGGGTG",fastq, max.mismatch=1,with.indels=FALSE))
  start_RPI<-as.vector(unlist(start_RPI))
  start<-append(start, start_RPI)
}
min<-(div-1)*200000-199999
max<-l
fastq<-seq[min:max]
start_RPI<-startIndex(vmatchPattern("TGGAATTCTCGGGTG",fastq, max.mismatch=1,with.indels=FALSE))
start_RPI<-as.vector(unlist(start_RPI))
start<-append(start, start_RPI) #this is the final vector with the start coordinate of RPI adapter (the border btw N9 and RPI)
rm("fastq")
length(start)

reads_trimmed<-narrow(x=seq, start = 1, end = (start-1)) 
plot(table(width(reads_trimmed)))



########################################################
########################################## 3.  Getting rid of adapter dimers and reads with too short inserts

reads_F<-reads_trimmed[width(reads_trimmed)>=40]
plot(table(width(reads_F)))
length(reads_F)

writeFasta(reads_F, "Forward_trimmed2.fasta") #write sequences with an insert into a Fasta file (each read starts with 4nt barcode followed by N11, insert and N9)









rm(list=setdiff(ls(), c("folder_sample","Barcode_4nt")))








############################   2*150
gc()
setwd(folder_sample)
setwd("2_150")
F_150<-as.character("KD403_FragSeq_Rep3_HTS2_2x150_Forward_reads") #insert the name of the forward reads
R_150<-as.character("KD403_FragSeq_Rep3_HTS2_2x150_Reverse_reads") #insert the name of the reverse reads


###forward reads
reads_F=readFastq(F_150)
length(reads_F) #how many reads
width_F<-as.numeric(width(reads_F))
stat_F<-as.data.frame(table(width_F))
rm(reads_F)


##reverse reads
gc()
reads_R=readFastq(R_150)
width_R<-as.numeric(width(reads_R))
stat_R<-as.data.frame(table(width_R))
rm(reads_R)


##we will analyze read quality by groups of reads with the same read length. Therefore to keep forward and reverse reads in the same order we will trimm each pair of forward and reverse reads so that each pair has the same read length
####min width for each pair
min_width<-pmin(width_F, width_R)
rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt", "min_width", "F_150", "R_150")))
gc()

##forward reads
reads_F=readFastq(F_150)
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
rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt", "min_width", "For", "R_150")))
writeFasta(For, "Forward_trimmed.fasta")  


#reverse reads
gc()
rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt", "min_width", "R_150")))


reads_R<-readFastq(R_150)
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
writeFasta(DNAStringSet(Combined_insert), "Combined_2_150.fasta")


rm(list=setdiff(ls(),c("folder_sample","Barcode_4nt")))




#######combine_all_together (1x150 and 2x150)
setwd(folder_sample)
setwd("1_150")


reads1<-append(readFasta("Forward_trimmed1.fasta"), readFasta("Forward_trimmed2.fasta"))
length(reads1)

setwd(folder_sample)
setwd("2_150")
reads2<-readFasta("Combined_2_150.fasta")




reads<-append(reads1, reads2)
length(reads) #number of all reads with an isert >=16


plot(table(width(reads)))
setwd(folder_sample)
writeFasta(reads, "Combined_reads")






###which contain 4nt barcode?
reads<-DNAStringSet(sread(reads))
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


###Here we will divide reads with identical sequence of the insert and molecular identifiers (N11 and N9) into groups and reduce amount of reads per group to just one by finding consensus sequence for each group
###This step will take several hours to be completed but you'll see the number of a loop just completed

reads<-table(reads)
reads<-data.frame(reads)
reads<-as.character(reads$reads)
length(reads)

plot(table(width(reads)))


width<-data.frame(table(width(reads)))

width_num<-as.numeric(as.character(width$Var1))
reads<-DNAStringSet(reads)
length(width_num) #how many loops


list<-list()
i=1
for (s in 1:length(width_num)){
  w<-width_num[s]
  r<-reads[width(reads)==w]
  
  while (length(r)>0) {
    pat<-as.character(r[1])
    hits<-vcountPattern(pat,r,
                        max.mismatch=0,with.indels=FALSE, fixed=FALSE)
    mult_reads<-r[hits==1]
    r<-r[hits!=1]
    list[[i]]<-mult_reads
    i=i+1
    
  }
  print(as.character(s))
}

length<-as.numeric(length(list))
length


ls<-lengths(list)
ls_1<-ls==1
list_short<-list[ls_1]


list_long<-list[!ls_1]
length(list_long)

norm_reads<- DNAStringSet(do.call(c,unlist(list_short)))
norm_reads<-as.character(norm_reads)

length_list<-length(list_long)
length_list

for (i in 1:length(list_long)) {
  inst<-list_long[[i]]
  seq<-consensusString(inst)
  norm_reads<-append(norm_reads, seq)
  print(paste(as.character(i), "/", sep=""))
}



reads<-DNAStringSet(norm_reads) #these are the reads with removed overamplified sequences (number of sequences with the same insert and N11+N9 is substituted with one read having a consensus sequence of this group of reads)

reads1<-substr(reads, 12, width(reads)-9) #removing N11 and N9
writeFasta(DNAStringSet(reads1), "reads_for_alignment.fasta") #these are the sequences of the inserts ready for alignment

plot(table(width(reads1)), main="Read length") #plot of read lengths

########################################################################
#########################################align reads to the genome

rm(list=setdiff(ls(), c("folder_sample")))


setwd(folder_sample) #make sure that your folder contains "kd403-real.fasta" reference genome 
insert<-readFasta("reads_for_alignment.fasta")
insert<-sread(insert)



#subject for mapping
mychr <- sread(readFasta("kd403-real.fasta")) # Creates sample chromosome
mychr_rev<-reverseComplement(mychr) #reverseComplement of the reference genome




insert<-gsub("R","N",insert)
insert<-gsub("Y","N",insert)
insert<-gsub("M","N",insert)
insert<-gsub("K","N",insert)
insert<-gsub("S","N", insert)
insert<-gsub("W","N",insert)
insert<-gsub("H","N", insert)
insert<-gsub("B","N", insert)
insert<-gsub("V","N", insert)
insert<-gsub("D","N", insert)
insert<-gsub("N", "a", insert)
insert<-DNAStringSet(insert)






width<-width(insert)
tab<-data.frame(table(width))

###mapping reads to KD403 genome
MM<-2 #mismatches allowed

unique_reads<-c()


for (i in 1:nrow(tab)){
  w<-as.numeric(as.character((tab[i,1])))
  
  set<-insert[width(insert)==w]
  mypdict <- PDict(set, max.mismatch=MM)
  
  F_amount <- countPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM)
  F_amount_rev <- countPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM)
  log1<- F_amount == "1" & F_amount_rev == "0"
  log2 <- F_amount == "0" & F_amount_rev == "1" 
  log11<-(log1==TRUE)|(log2==TRUE)
  table(log11)
  
  set_unique<-set[log11]
  unique_reads<-append(unique_reads, set_unique)
}


length(unique_reads) #amount of reads uniquely aligned to KD403 genome


##Determine position of alignment for each uniquely aligned read
insert_un<-unique_reads
width<-width(insert_un)
tab<-data.frame(table(width))

w<-as.numeric(as.character((tab[1,1])))
set<-insert_un[width(insert_un)==w]
mypdict <- PDict(set, max.mismatch=MM)

F_amount <- countPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM)
F_amount_rev <- countPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM)
mysearch <- matchPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
mysearch_rev<-matchPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM) # Searches all dictionary entries against chromosome.

Forward<-data.frame(mysearch)
Reverse<-data.frame(mysearch_rev)

log<-F_amount==1
reads<-set[log]
Forward<-cbind(Forward, reads)

log<-F_amount_rev==1
reads<-set[log]
Reverse<-cbind(Reverse, reads)


for (i in 2:nrow(tab)){
  w<-as.numeric(as.character((tab[i,1])))
  
  
  set<-insert_un[width(insert_un)==w]
  mypdict <- PDict(set, max.mismatch=MM)
  
  F_amount <- countPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM)
  F_amount_rev <- countPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM)
  
  mysearch <- matchPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
  mysearch_rev<-matchPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
  
  Forward1<-data.frame(mysearch)
  Reverse1<-data.frame(mysearch_rev)
  
  
  log<-F_amount==1
  reads<-set[log]
  Forward1<-cbind(Forward1, reads)
  
  log<-F_amount_rev==1
  reads<-set[log]
  Reverse1<-cbind(Reverse1, reads)
  
  
  Forward<-rbind(Forward, Forward1)
  Reverse<-rbind(Reverse, Reverse1)
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

