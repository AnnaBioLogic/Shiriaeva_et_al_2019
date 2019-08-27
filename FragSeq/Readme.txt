This is a custom code for FragSeq analysis of short fragments purified from E.coli cells with type I-E or type I-F CRISPR-Cas system.
Some of the fastq files are large. I analyzed all data on my laptop with 16GB of RAM and it was enough.
The analysis is performed in R, ShortRead package will be required.
For some samples we repeated sequencing of the same NGS libraries several times so there are samples like KD518_rep3 for which 1x150, 2x75 and 2x150 sequencing was performed.
Make a separate folder for each sample on your computer (for example, "Folder X") and inside of that folder make folders "1_150", "2_75" or "2x150" for each type of sequencing (these is how it is organised on my computer and how it appears in the code).
Place unpacked fastq files into corresponding "1_150", "2_75" or "2x150" folders. Processed reads will be saved in the folder "Folder X"
Place the reference sequence fasta files into the "Folder X". After alignment is complete, the alignment table will be saved in the same folder.
