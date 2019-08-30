###Install these libraries:

source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
biocLite("Rsamtools")
biocLite("ggplot2")


rm(list=ls())
library(ShortRead)




#setwd to a folder that contains  .fastq files; 
setwd("E:/Genomic DNA/4.3")






reads_F<-readFastq("Strain_comparison_experiment_KD518_WithCasInduction_GenomicDNA_NEBNext_UltraII_Rep3.fastq") #reads_F contains all reads
length(reads_F)
boxplot(width(reads_F))
table(width(reads_F))

width_F<-as.numeric(width(reads_F))
stat<-as.data.frame(table(width_F))
stat$Freq<-as.numeric(as.character(stat$Freq))
stat$width_F<-as.numeric(as.character(stat$width_F))



PQcutoff <- 20 #we will substitute bases with PhredQuality less than this to N
rm(width_F)



M1_R1<-c()
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
    M1_R1<-append(M1_R1, seq_R)  
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
        M1_R1<-append(M1_R1, seq_R)        }
      
      
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
        M1_R1<-append(M1_R1, seq_R)
      }
      
    }
    
  }
  
}

rm(list=setdiff(ls(), "M1_R1"))

##how many Ns per read? we will keep no more than 2 Ns per read
hits_N=vcountPattern("N",M1_R1,
                          max.mismatch=0,with.indels=FALSE)
table(hits_N)
log<-hits_N<=2

good_reads<-M1_R1[log]
bad_reads<-M1_R1[!log]

start_N<-startIndex(vmatchPattern("N",bad_reads,max.mismatch=0, with.indels=FALSE))
first_N<-lapply(start_N, function(x) x[3])
first_N<-as.vector(unlist(first_N))

bad_reads_trimmed<-narrow(bad_reads, 1, first_N-1)

boxplot(width(bad_reads), main="Bad_reads")
boxplot(width(bad_reads_trimmed),main="Bad_reads_trimmed" )

reads<-c(good_reads, bad_reads_trimmed)

###cut 3' adapter from reads
hits_3prime=vcountPattern("AGATCGGAAGAGCA",reads,
                          max.mismatch=1,with.indels=FALSE)
log<-hits_3prime==1
table(log)


M1_R1_1<-reads[log]
M1_R1_2<-reads[!log]

if(length(M1_R1_1)>=1){
start_adapter<-startIndex(vmatchPattern("AGATCGGAAGAGCA",M1_R1_1,max.mismatch=1, with.indels=FALSE))
start_adapter<-as.vector(unlist(start_adapter))

M1_R1_1<-narrow(M1_R1_1, 1, start_adapter-1)
reads<-append(M1_R1_1, M1_R1_2)

}

Rpattern <- "AGATCGGAAGAGCA"
reads<- trimLRPatterns( Rpattern = Rpattern, subject = reads)

log<-width(reads)>=20 #keep the readslonger than 20 nt
reads<-reads[log]
rm(list=setdiff(ls(), "reads"))

dir.create("processed_reads")
setwd("processed_reads")
writeFasta(reads ,file="reads.fasta")

##reads.fasta is further used for alignment to the reference genome (file named "kd403-real.fasta") in UniproUgene program http://ugene.net
#In Unipro Ugene use Bowtie2 to map the reads (end-to-end alignment; 1 mismatch allowed) to "kd403-real.fasta genome"
# The program will write down the file called "kd403-real.sam.bam_sorted.bam" into the same folder where your reads are. Once you have it, go back to this R code for further analysis








############################################################################
#######################################################
##################################Extracting BAM files with mapq==42
library(Rsamtools)



rm(list=ls())
bamFile<-BamFile("kd403-real.sam.bam_sorted.bam")
## filter to a single file
####filter all alignments with mapq==42 to bam_42_all.bam file
p3<-ScanBamParam(what=scanBamWhat(), mapqFilter=42)
filter <- FilterRules(list(mapq = function(x) (x$mapq) == 42))
destination<-"bam_42_all.bam"
dest <- filterBam(bamFile, destination, param=p3, filter=filter)
countBam(dest)  



mychr <- sread(readFasta("kd403-real.fasta")) # Creates sample chromosome.
width<-width(mychr)



####alignment normalization

alignsTest1<-readGAlignments('bam_42_all.bam')
length(alignsTest1) #(this is is how many reads aligned you have)
COV1<-coverage(alignsTest1)
UNCOV1<-unlist(COV1$KD403)
UNCOV1_L<-runLength(UNCOV1)
UNCOV1_V<-runValue(UNCOV1)
UNCOV1_vector<-rep(UNCOV1_V,UNCOV1_L )

name="KD518_rep.3"


##Whole genome plot
plot(UNCOV1_vector, type="h", xlim=c(1, width(mychr)), ylim=c(0, 500), xlab="bp", ylab="Coverage", main=name )

#################plot PPS-/+300kb

plot(UNCOV1_vector, type="h", xlim=c(275539, 875570), ylim=c(0, 500), xlab="bp", ylab="Coverage", main=name )
lines(cbind(575540,-10000),col="red", type="h", lwd=1.5) 
lines(cbind(575540,10000),col="red", type="h", lwd=1.5) 
lines(cbind(575570,-10000),col="red", type="h", lwd=1.5) 
lines(cbind(575570,10000),col="red", type="h", lwd=1.5) 





#####calculate mean coverage in different genomic regions;


###how much in -200kb....PPS...+100kb

min<-575538-200000  #coordinate of 200kb upstream of the PPS
max<-575570+100000 #coordinate of  100kb downstream of the PPS
middle<-575538+(575570-575538)/2 #coordinate of the PPS

length(UNCOV1_vector)
width(mychr)


mean(UNCOV1_vector)  #mean coverage - whole genome


min200<-mean(UNCOV1_vector[(min-5000):(min+5000)])
max100<-mean(UNCOV1_vector[(max-5000):(max+5000)])
flanks<-(min200+max100)/2




PPS<-mean(UNCOV1_vector[(middle-5000):(middle+5000)])






terC<-3027258
ter<-mean(UNCOV1_vector[(terC-5000):(terC+5000)])

###You can find following values in Supplementary Table 3:
flanks  # mean coverage - PPS-flanking regions
PPS #mean coverage - PPS
ter # mean coverage - terC
flanks/PPS
flanks/ter









##################################################################################
#Here we will remove repeated regions from the alignment and calculate normalized coverage per each 1 kb of genome
#make sure that there are files "repeats_of_KD403" and "kd403-real.fasta" in the folder 

mychr <- sread(readFasta("kd403-real.fasta")) # Creates sample chromosome.
width<-width(mychr)


repeats<-read.table("repeats_of_KD403", row.names = NULL)



my_table<-data.frame(c(1:width), UNCOV1_vector)
###remove repeated elements


for (i in 1:nrow(repeats)){
  start<-as.numeric(as.character(repeats[i,2]))
  end<-as.numeric(as.character(repeats[i,3]))
  
  
  log<-my_table$c.1.width.>=start&my_table$c.1.width.<=end
  my_table[log, 2]=0
  
  
}



####coverage in  bins

non_zero<-(my_table$UNCOV1_vector != 0)
my_table_non_zer0<-data.frame(my_table$c.1.width.[non_zero], my_table$UNCOV1_vector[non_zero])
colnames(my_table_non_zer0)=c("Coordinate", "Coverage")

my_table<-my_table_non_zer0



w<-1000 ####bin
###E coli genome 4624855

k<-width/w
k<-ceiling(k)
avcov<-c()
x<-c()


for (i in 1:k) {
  x1<-i*w-w
  y1<-i*w
  
  
  
  my_table_local<-my_table[((as.numeric(my_table[,1])<=y1)&(as.numeric(my_table[,1])>x1)),]
  avcov1<-(sum(as.numeric(as.character(my_table_local[,2])))/nrow(my_table_local))
  
  avcov<-c(avcov, avcov1)
  coor<-(((y1-x1)/2)+x1)
  x<-c(x, coor)
}


df<-data.frame(x, avcov)
df <- replace(df, is.na(df), 0)
log<-df$avcov!=0
df<-df[log,]
df <- setNames(df, c("middle coordinate","average coverage" ))
sum<-sum(as.numeric(as.character(df$`average coverage`)))
coverage_normalized_percent<-100*(as.numeric(as.character(df$`average coverage`)))/sum
df<-cbind(df, coverage_normalized_percent)


write.table(df,file='genomic_cov_per_1kb')


plot(df[,1], df[,2])
library(ggplot2)

df2<-df[df$`middle coordinate`>=1 & df$`middle coordinate`<1000000,]
ggplot(df, aes(df$`middle coordinate`, df$coverage_normalized_percent))+geom_point()+geom_smooth(span=0.2)
ggplot(df2, aes(df2$`middle coordinate`, df2$coverage_normalized_percent))+geom_point()+geom_smooth(span=0.15)



