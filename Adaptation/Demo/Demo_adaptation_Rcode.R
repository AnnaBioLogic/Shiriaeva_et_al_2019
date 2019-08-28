source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")


library(ShortRead)

rm(list=ls())
setwd("D:/Adaptation/Demo") #specify folder with fastq files; make sure that the folder also contains reference genome sequence "kd403-real.fasta"

ID="Demo_adaptation"

reads_F<-readFastq("Adaptation _demo.fastq")
length(reads_F)
  width_F<-width(reads_F)
  
stat<-as.data.frame(table(width_F))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$width_F<-as.numeric(as.character(stat$width_F))



PQcutoff <- 20 #we will substitute bases with PhredQuality less than this to N


For<-c()
for (i in 1:nrow(stat)){
  width<-stat[i,1]
  reads<-reads_F[width(reads_F)==width]
  
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


M1<-For

##How many reads contain CRISPR>=1

length(M1[vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1,max.mismatch=2,with.indels=F)>=1])+
  length(M1[vcountPattern("GGTTTATCCCCGCTGGCGCGGGGAACTC",M1,max.mismatch=2,with.indels=F)>=1])

###How many reads contain CRISPR>=2
length(M1[vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1,max.mismatch=2,with.indels=F)>=2])+
  length(M1[vcountPattern("GGTTTATCCCCGCTGGCGCGGGGAACTC",M1,max.mismatch=2,with.indels=F)>=2])


###How many reads contain CRISPR==3

length(M1[vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1,max.mismatch=2,with.indels=F)==3])+
  length(M1[vcountPattern("GGTTTATCCCCGCTGGCGCGGGGAACTC",M1,max.mismatch=2,with.indels=F)==3])



M1_dir<-M1[!(vcountPattern("GGTTTATCCCCGCTGGCGCGGGGAACTC",M1,max.mismatch=2,with.indels=F)>=1)&vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1,max.mismatch=2,with.indels=F)>=2]
length(M1_dir)

M1_rev<-M1[vcountPattern("GGTTTATCCCCGCTGGCGCGGGGAACTC",M1,max.mismatch=2,with.indels=F)>=2&!(vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1,max.mismatch=2,with.indels=F)>=1)]
length(M1_rev)


M1_rev1<-chartr("ATGC", "TACG", M1_rev)# Returns complement for given DNA sequences.
x<-strsplit(as.vector(M1_rev1), "") # Vectorizes many sequences.
x<-lapply(x,rev)#reverses vectore
M1_rev1<-sapply(x,paste, collapse="")#collapse vectors to strings.


M1_all<-append(M1_rev1,as.vector(M1_dir))
M1_all<-DNAStringSet(M1_all)

length(M1_all)

hist(width(M1_all))#### 

##arrays with 2 sps
M1_2<-M1_all[vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1_all,max.mismatch=2,with.indels=F)==2]

##arrays with 3 sps
M1_3<-M1_all[vcountPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1_all,max.mismatch=2,with.indels=F)==3]

##arrays with 2 sps
start_2<-startIndex(vmatchPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1_2,max.mismatch=2,with.indels=F))

end_2<-endIndex(vmatchPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1_2,max.mismatch=2,with.indels=F))

start_2<-as.vector(unlist(start_2))

start_2<-start_2[seq(0,length(start_2),2)]

end_2<-as.vector(unlist(end_2))

end_2<-end_2[seq(1,length(end_2),2)]

spacers_2<-substr(M1_2,end_2+1,start_2-1)

###########################
#########################
######arrays with three sp

start_3<-startIndex(vmatchPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1_3,max.mismatch=2,with.indels=F))

end_3<-endIndex(vmatchPattern("GAGTTCCCCGCGCCAGCGGGGATAAACC",M1_3,max.mismatch=2,with.indels=F))

start_3<-matrix(as.vector(unlist(start_3)),nrow=length(M1_3), ncol=3,byrow=TRUE)


end_3<-matrix(as.vector(unlist(end_3)),nrow=length(M1_3), ncol=3,byrow=TRUE)

spacers_31<-substr(M1_3,end_3[,1]+1,start_3[,2]-1)
spacers_32<-substr(M1_3,end_3[,2]+1,start_3[,3]-1)


spacers_all<-spacers_2
spacers_all<-append(spacers_all, spacers_31)
spacers_all<-append(spacers_all, spacers_32)

length(spacers_all) #how many spacers

hist(width(spacers_all),breaks=100)

###remove a preexisting spacer 
old<-"TCAAACAACCGACCTTGTTGTTCGCTATT"
spacers_all<-DNAStringSet(spacers_all)
hits=vcountPattern(old,spacers_all,
                       max.mismatch=3,with.indels=FALSE)

log<-hits==1
table(log)
spacers_all<-spacers_all[!log]
length(spacers_all) #new spacers without the preexisting spacer
table(width(spacers_all))

sp_33<-spacers_all[width(spacers_all)==33] #these are extracted 33 bp spacers (we will align them to "kd403-real.fasta" genome)
sp_32<-spacers_all[width(spacers_all)==32]#these are extracted 32 bp spacers
sp_34<-spacers_all[width(spacers_all)==34]#these are extracted 34 bp spacers

df<-data.frame(table(width(spacers_all)))

write.table(df, "all_spacers_lengths")


writeFasta(DNAStringSet(sp_33),file="spacers_33.fasta")
writeFasta(DNAStringSet(sp_32),file="spacers_32.fasta")
writeFasta(DNAStringSet(sp_34),file="spacers_34.fasta")









new_spacers<-sread(readFasta("spacers_33.fasta"))
length(new_spacers) #how many spacers of 33 bp length

MM<-2 #mismatches allowed



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
new_spacers<-gsub("N", "a", new_spacers)
new_spacers<-new_spacers[width(new_spacers)==33]
length(new_spacers)



#subject for mapping
mychr <- sread(readFasta("kd403-real.fasta")) # Creates sample chromosome.
#creates rev sample chromosome
mychr_rev<-chartr("ATGC", "TACG", mychr)# Returns complement for given DNA sequences.
x<-strsplit(as.vector(mychr_rev),"") # Vectorizes many sequences.
x<-lapply(x,rev)#reverses vectore
mychr_rev<-sapply(x,paste, collapse="")#collapse vectors to strings.

#counting number of alignment for each spacer

mypdict <- PDict(new_spacers[width(new_spacers)==33], max.mismatch=MM)
amount <- countPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM)
amount_rev <- countPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM)
n <-c(new_spacers)
df<-data.frame(n, amount, amount_rev)
df$sum=df$amount+df$amount_rev
table(df$sum>0) #how many spacers were aligned at least once

#select only spacers with unique alignment
df1 <- subset(df, amount == "1" & amount_rev == "0")
df2 <- subset(df, amount == "0" & amount_rev == "1")
df3 <- rbind(df1, df2)
spacers <- df3[,1]
new_spacers <- as.character(spacers)
length(new_spacers) #how many spacers were uniquely aligned



#we will map only uniquely aligned spacers

mypdict <- PDict(new_spacers[width(new_spacers)==33], max.mismatch=MM)
mysearch <- matchPDict(mypdict, DNAString(as.character(mychr)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
mysearch_rev<-matchPDict(mypdict, DNAString(as.character(mychr_rev)), max.mismatch=MM) # Searches all dictionary entries against chromosome.
sum(elementNROWS(mysearch))#how many on dir strand
sum(elementNROWS(mysearch_rev))#how many on rev strand

############making a table with coordinates of read alignment and PAM sequence
P_start<-(-2) #Indicate PAM start position with respect to the spacer start (If=(-2), it means that PAM starts 2 letters upstream of the spacer start coordinate)
P_end<-(-0) #Indicate PAM end position with respect to the spacer start (If=(0), it means that PAM end position is the same as the spacer start coordinate like in E.coli type I-E)
spacer_length<-33



dir_11<-as.matrix(table(unlist(startIndex(mysearch))))#start coordinate on dir strand

PAM_start<-as.numeric(rownames(dir_11))+P_start
PAM_end<-as.numeric(rownames(dir_11))+P_end
spacer_start<- as.numeric(rownames(dir_11))
spacer_end<- as.numeric(rownames(dir_11))+spacer_length-1

a2<-c(2)
for (i in seq(1, nrow(dir_11)))
{
  a2[i]<-substr(mychr,PAM_start[i],PAM_end[i])
}



a3<-c(2)
for (i in seq(1, nrow(dir_11)))
{
  a3[i]<-substr(mychr,spacer_start[i],spacer_end[i])
}

dir_11<-cbind(dir_11,a2)  # 
dir_11<-cbind(dir_11,a3)# 
dir_11<-cbind(dir_11, "1")
colnames(dir_11)<-c("quantity", "PAM","spacer", "strand")


###############
############the same for the rev strand
rev_11<-as.matrix(table(unlist(startIndex(mysearch_rev))))#start position


PAM_start<-as.numeric(rownames(rev_11))+P_start
PAM_end<-as.numeric(rownames(rev_11))+P_end
spacer_start<- as.numeric(rownames(rev_11))
spacer_end<- as.numeric(rownames(rev_11))+spacer_length-1

a2<-c(2)
for (i in seq(1, nrow(rev_11)))
{
  a2[i]<-substr(mychr_rev,PAM_start[i],PAM_end[i])
}



a3<-c(2)
for (i in seq(1, nrow(rev_11)))
{
  a3[i]<-substr(mychr_rev,spacer_start[i],spacer_end[i])
}

rev_11<-cbind(rev_11,a2)  
rev_11<-cbind(rev_11,a3)
rev_11<-cbind(rev_11, "2")
colnames(rev_11)<-c("quantity", "PAM","spacer", "strand")
rownames(rev_11)<-width(mychr)-as.numeric(rownames(rev_11))+1



dir_rev_11<-rbind(dir_11, rev_11) #this is a combined table for both strands ("1"=direct strand; "2"=reverse strand)


write.table(dir_rev_11,file='sp33_aligned_to_genome') #the table includes information on spacer sequence, position of the first nucleotide, PAM and quantity of the spacer



###############Pictures

my_table<-read.table('sp33_aligned_to_genome', row.names=NULL)
colnames(my_table)<-c("position", "quantity", "PAM", "spacer", "strand")
my_table$position<-as.numeric(as.character(my_table$position))


sum<-sum(as.numeric(as.character(my_table$quantity)))
percent<-100*(as.numeric(as.character(my_table$quantity)))/sum
my_table<-cbind(my_table, percent)

dir<-my_table[my_table$strand=="1",]
rev<-my_table[my_table$strand=="2",]
rev[,2]<-(-1)*as.numeric(as.character(rev[,2]))
rev[,6]<-(-1)*as.numeric(as.character(rev[,6]))


###WHOLE GENOME 
plot(dir[,c(1, 6)], type="h", xlim=c(0, width(mychr)), ylim=c(min(rev[,6]), max(dir[,6])),  ylab="Percent of spacers", xlab="bp",   cex.lab=1.2, main=ID)
lines(rep(0,width(mychr),col=1))
lines(dir[,c(1,6)],col="blue",type = "h", lwd=1.7)
lines(rev[,c(1,6)],col="red",type = "h", lwd=1.7)



lines(cbind(575540,-40000),col="black", type="h", lwd=1.5) #PPS
lines(cbind(575540,40000),col="black", type="h", lwd=1.5) #PPS
lines(cbind(575570,-40000),col="black", type="h", lwd=1.5) #PPS
lines(cbind(575570,40000),col="black", type="h", lwd=1.5) #PPS


