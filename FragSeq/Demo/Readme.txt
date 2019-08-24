To demo this code we will use a subset of data for KD518_FragSeq_Rep1 experiment.
The whole analysis is performed in RStudio, ShortRead package will be required.
Make a folder for demo testing ("Demo") that contains a folder called "2_75" with fastq files: "KD518_FragSeq_Rep1_HTS2_2x75_Forward_reads2" and "KD518_FragSeq_Rep1_HTS2_2x75_Reverse_reads2" that contain 3985332 reads each.
Place the reference genome sequence "kd403-real.fasta" into the "Demo" folder.


Run Lines 8-203: reads are trimmed and N is placed at the positions of low PhredQuality score. Processed forward and reversed reads are written into your 2x75 folder in fasta format.
Lines 209-213: remove all files except for the pathway to the folder with your reads ("folder_sample") and the 4nt barcode sequence "Barcode_4nt" that we specified in the very beginning (in our HTS run- different strain were marked with different 4nt barcodes that are first 4 nt in the forward reads).
Upload Forward_trimmed.fasta and Reverse_trimmed.fasta files, reverseComplement returns reverse complement sequence of Reverse_trimmed.fasta reads. So now if For and Rev reads are overlapping, a part of their sequences is identical.

Lines 216-233: testing if For and Rev reads overlap. If so, we will get combined sequence that starts with the first nucleotide of the For read and ends with the nucleotide complementary to the first nucleotide of the Rev read.
All combined reads can be found in "Combined_insert" 
Line 236: plot the lengths of combined reads.

The adapters used for FragSeq were designed in a way that reads in "Combined_insert" file start with 4nt barcode followed by 11nt unique molecular identifier (N11 UMI) followed by insert (sequence of the DNA fragment) follwed by 9nt unique molecular identifier (N9 UMI).
Because of the small size of the insert it is hard to get rid of adapter dimers so we have plenty of them in our libraries.
Lines 240-243: remove adapter dimers and all reads with insert shorter than 16nt. 

Lines 265-275: keep only those reads that start with right 4nt barcode. Than remove this 4nt barcode from your reads.
For consistency between all samples we will analyze fragments 16-100 nt in all strains.

If an insert is abundant in your sample you will have many reads with the same insert but different N11 and N9 UMIs. If the same combination of an insert, N9 and N11 is sequenced many times most likely these reads are PCR duplicates. We will remove all but one read of identical sequences.
The read that we will keep for each group of identical N11_insert_N9 sequences is a consensus sequence of all sequences in this group (Lnes 285-352).
Than we will remove N11 and N9 UMIs and write fragment sequences into a file "reads_for_alignment.fasta".
Line 354: plotting lengths of the fragment sequences ready for alignment

Lines 362-549: mapping inserts to "kd403-real.fasta" reference genome. Only reads with unique alignment are analyzed.
Table "Alignment_table_both_strands" will appear in your "Demo" folder. This table includes information about start and end coordinate of each read (for reads derived from the reverse strand of KD403, the start coordinate in this table is actually the coordinate of fragments' 3' ends).
The table also contains information about fragment lengths, strand where fragments were derived from and fragment sequences.
Along with the sequences of fragments there are sequences of each fragment plus 20 nt upstream of fragment's 5' end in the chromosome and 20 nt downstream of fragment's 3' end (the column called ""(20nt_up_chromosome)+READ+(20nt_dw_chromosome)").

Overall, the demo testing time for "KD518_FragSeq_Rep1_HTS2_2x75_Forward_reads2" and "KD518_FragSeq_Rep1_HTS2_2x75_Reverse_reads2" on my laptop is about 20 minutes.
17934 reads were uniquely aligned to the genome.
