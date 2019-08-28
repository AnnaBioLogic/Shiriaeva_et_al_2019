For demo version use "Adaptation _demo.fastq".
Make a folder containing the "Adaptation _demo.fastq" file and a reference genome sequence "kd403-real.fasta" 
The "Adaptation _demo.fastq" file has 27383 reads; 19988 reads have one CRISPR repeat; 9346 reads have two CRISPR repeats; 0 reads have 3 CRISPR repeats.
After removal of a preexisting self-targeting spacer you will get 9345 newly acquired spacer sequences.
9188 new spacers have a length of 33bp. We will align these spacers to the "kd403-real.fasta" genome sequence (with two mismatches allowed).
8360 are aligned to the genome. 8197 spacers are uniquely aligned. We will further analyze only uniquely aligned spacers.
A table named "sp33_aligned_to_genome" will be created in the same folder where your fastq file is. 
The table contains information about the position of the first protospacer nucleotide (last nucleotide in the case of spacers mapped to the reverse strand); spacer sequence; its quantity; PAM sequence and strand (1 - direct strand; 2 - reverse strand).
The "sp33_aligned_to_genome" for the "Adaptation _demo.fastq" has 1702 rows. 
You will get a picture where X axis corresponds to E.coli KD403 genome; Y axis - percent of protospacers. Above 0 - spacers mapped to the top nontarget strand (blue); below - spacers mapped to the bottom target strand (red).
Run time for "Demo_adaptation_Rcode.R" is about 1 minute. For larger fastq files like "Primed_adaptation_KD403_rep1" which is 1.8 Gb it takes ~20 minutes to complete the code on my laptop.
