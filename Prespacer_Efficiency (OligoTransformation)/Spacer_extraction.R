source("https://bioconductor.org/biocLite.R") 
biocLite("ShortRead")#install the package




library(ShortRead) #load the package


rm(list=ls())

setwd("D:/OneDrive/Results/Oligo transformation/HTS_new/1_2") #set way to the folder with a fastq file
reads=readFastq("G33_C33_OligoTransformation_Rep1") #upload the fastq file
length(reads) ###initial number of reads

dir.create("PQ14 0N_longer105")
setwd("D:/OneDrive/Results/Oligo transformation/HTS_new/1_2/PQ14 0N_longer105")

width<-as.numeric(width(reads))
barplot(table(width)) #shows distribution of read lengths
log<-width>=105 #check if the read length>=105nt (Leader(54nt)+Repeat(28nt)+1st spacer (at least 23nt))
table(log)#TRUE = how many reads have length>=105nt
reads_long<-reads[log] #keep only reads >=105nt



###Determine base quality at each position, subsitute nucleotides with quality less than the quality cutoff (PQcutoff) with "N" (We used PQcutoff=14 for this HTS run) 
PQcutoff<-14
length(reads_long)
width_F<-width(reads_long)

stat<-as.data.frame(table(width_F))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$width_F<-as.numeric(as.character(stat$width_F))

For<-c()
for (i in 1:nrow(stat)){
  width<-stat[i,1]
  reads<-reads_long[width(reads_long)==width]
  
  l<-length(reads)
  
  if (l<=300000) {
    reads_F1<-reads
    seq_R1 <- sread(reads_F1) 
    q <- PhredQuality(quality(quality(reads_F1))) 
    MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
    
    positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
    N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
    N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet") 
    rm(reads_F1)
    gc()
    seq_R <- replaceLetterAt(seq_R1, positions, N_tr) 
    For<-append(For, seq_R)  
  }
  else {
    whole<-l%/%300000
    rest<-l%%300000
    
    for (b in 1:(whole+1)) {
      if (b<(whole+1)) {
        reads_F1<-reads[(b*300000-299999):(b*300000)]
        seq_R1 <- sread(reads_F1) 
        q <- PhredQuality(quality(quality(reads_F1))) 
        MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
        
        positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
        N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
        N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet") 
        rm(reads_F1)
        gc()
        seq_R <- replaceLetterAt(seq_R1, positions, N_tr) 
        For<-append(For, seq_R)}
      
      
      else {
        reads_F1<-reads[(b*300000-299999):length(reads)]
        seq_R1 <- sread(reads_F1) 
        q <- PhredQuality(quality(quality(reads_F1))) 
        MATR <- matrix(charToRaw(as.character(unlist(q))), nrow=length(q), byrow=TRUE) 
        
        positions <- MATR < charToRaw(as.character(PhredQuality(as.integer(PQcutoff)))) 
        N <- DNAString(paste(rep.int("N", width(seq_R1)[1]), collapse="")) 
        N_tr <- as(Views(N, start=1, end=rowSums(positions)), "DNAStringSet") 
        rm(reads_F1)
        gc()
        seq_R <- replaceLetterAt(seq_R1, positions, N_tr) 
        For<-append(For, seq_R)
      }
      
    }
    
  }
  
}
M1<-For #M1 contains your reads

rm(list=setdiff(ls(),c("M1")))


M1<-gsub("R","N",M1)
M1<-gsub("Y","N",M1)
M1<-gsub("M","N",M1)
M1<-gsub("K","N",M1)
M1<-gsub("S","N",M1)
M1<-gsub("W","N",M1)
M1<-gsub("H","N",M1)
M1<-gsub("B","N",M1)
M1<-gsub("V","N",M1)
M1<-gsub("D","N",M1)


table(vcountPattern("N",M1,max.mismatch=0,with.indels=F)) #shows how many reads have 0, 1, 2.... Ns per read
log<-vcountPattern("N",M1,max.mismatch=0,with.indels=F)<=0 #we will keep reads with no more than 0 Ns
table(log) #no more than 0N per read
M1<-DNAStringSet(M1[log])


#####How many reads contain at least 1 repeat?
rep="GTGTTCCCCGCGCCAGCGGGGATAAACC"
LOG1<-(vcountPattern(rep,M1,max.mismatch=3,with.indels=F)>0)
M1_r<-M1[LOG1] ###the reads with repeats
M1_NOT_r<-M1[!LOG1] ###the reads without repeats
length(M1_r)



#where does the first DR stop??
end<-endIndex(vmatchPattern(rep,M1_r,max.mismatch=3,with.indels=F))
end<-unlist(lapply(end, function(x) x[1])) #take the end coordinate of the first repeat only

log<-(end+1)<(width(M1_r))
M1_r<-M1_r[log]
end<-end[log]

M1_trimmed<-narrow(M1_r, end+1, width(M1_r))  ###trimm reads so they now start with the first nucleotide of a spacer (G in the case of a typial E.coli spacer)
log<-width(M1_trimmed)>=24#(at least 24 nucleotides of the first spacer)
table(log)
M1_trimmed<-M1_trimmed[log]


#####how many reads contain second repeat? (these are the reads with one additional spacer)
LOG1<-vcountPattern(rep,M1_trimmed,max.mismatch=3,with.indels=F)==1
table(LOG1)

adapted<-M1_trimmed[LOG1] ###these are the reads with additional repeat
nonadapted<-M1_trimmed[!LOG1]
nonadapted_24<-narrow(nonadapted, 1,24)


#####How many nonadapted reads contain preexisting spacer?
sp0<-"GAGCACAAATATCATCGCTCAAAC"
log_sp<-vcountPattern(sp0,nonadapted_24,max.mismatch=2,with.indels=F)==1
table(log_sp)
length(log_sp[log_sp==TRUE]) #Nonadapted (Contain old spacer)
length(log_sp[log_sp!=TRUE]) #Nonadapted (Smth else instead of the old spacer)

maybe_adapted<-nonadapted_24[!log_sp]
length(adapted) #these are reads with one additional spacer
length(maybe_adapted) #not clear what follows the repeat
length(nonadapted) #nonadapted


##determine the position of the second repeat and extract new spacers
start<-startIndex(vmatchPattern(rep,adapted,max.mismatch=3,with.indels=F))
start<-as.vector(unlist(start))
sp1<-narrow(adapted, 1, start-1)  #trimm by the start coordinate of the second repeat to keep just the spacer sequence
spacers<-sp1 #these are the sequences of new spacers



hist(width(spacers)) #shows distribution of spacer lengths
table(width(spacers))#shows distribution of spacer lengths

writeFasta(DNAStringSet(spacers),file="spacers.fasta") #writes spacer sequences into the file "spacers.fasta"
writeFasta(DNAStringSet(maybe_adapted),file="maybe_spacers.fasta") #writes spacer sequences into the file "maybe_spacers.fasta"

