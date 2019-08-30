For running the demo code use "Demo_total_genomic_DNA.fastq" file. The file has  1153489 reads.
Run the code to the line number 152 (takes less than a minute on my computer). A folder "processed_reads" will be created inside the folder with the fastq file. 
Place "kd403-real.fasta" reference sequence into the "processed_reads" folder.

Run Unipro UGENE. Go to "Tools"->"NGS data analysis" ->"Map reads to reference".
Choose "Bowtie2" as an alignment method. Specify a pathway to your reference sequence ("kd403-real.fasta").
Add "reads" file from the "processed_reads" folder.
Mode = end-to-end
Number of mismatches = 1
Library = Single-end
Press "Start"
Overall, it takes 1.5 min to map the Demo reads to the genome on my laptop. Once the procedure is complete close Ugene and go back to the R script.

Write the name of your bam file in line 173. Run lines 174-211 of the code (takes ~3 min for the Demo reads). Only reads with unique alignment will be kept.
You will get a picture of alignment near the PPS (red line=PPS). Not normalized coverage for each coordinate along the genome is shown.


Lines 217-255: calculate mean coverage in different genomic regions. For example, for the Demo sample mean coverage throughout the genome is 8.156947; mean coverage in the PPS-flanking regions is 8.118588 (The reads are from KD403 without cas inducers); mean coverage near the terC site is 7.080392. 

Run the rest of the code (~7 min). It will write a file called 'genomic_cov_per_1kb' with calculated normalized coverage per 1kb of the genome (mean coverage per each 1kb *100/sum of means).
Finally, you will get a scatter plot of normalized genomic coverage where each dot corresponds to the mean coverage per 1 kb with a smooth curve added. 
