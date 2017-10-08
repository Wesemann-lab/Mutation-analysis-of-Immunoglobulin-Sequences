# Mutation analysis of Immunoglobulin Sequences 

A tool for mutation analysis of high throughput monoclonal immunoglobulin sequences from Germinal center B-cells. 

This is the pipeline was used for analysis of data in the publication "Stochasticity enables BCR-independent germinal center initiation and antibody affinity maturation".

It is a very simple workflow which works mostly in R. Pear, Muscle and Jalview were used in intermediate steps. The functions used and the steps followed are descripted below.

Requirements:

1) R and R studio
2) Pear (Paired-End reAd mergeR) (https://sco.h-its.org/exelixis/web/software/pear/doc.html)
3) Muscle (http://www.drive5.com/muscle/)
4) Jalview (http://www.jalview.org/)


Make sure that all the bio-packages listed at the beginning of each R script are installed in R. 

Work flow:

 1) Merge the paired end reads by running Pear with the following parameters
		i)  minimum overlap for merge (-v) = 20
		ii) quality score threshold for trimming the low quality part of a read (Phred score) (-q) = 22
		iii)minimum possible length of the assembled sequences = 380 for heavy chain and 350 for light chain

 2) The following function in R replace intermittent low quality reads(Phred score <= 20) by “N” and converts Fastq file to a Fasta file. It also trims the sequences it finds four consecutive low quality reads. 
	Seq_File_filter.R - The input for this file is the output file from pear (merged reads) and a file containing the Phred scores along with their ascii characters(available with the files) to decode the characters into scores. This function also calls another function “Lower_score_N.R”. Both the functions should be sources into the workspace before running the first function.

 3) In-order to avoid very poor quality sequences, a filter is applied to allow sequences with a maximum of 2 percent ’N’s in them to be included in 	the next step. All sequences with greater than 2 % “N”s were excluded. 
	Percent_N_apply_file - The input for this function is the output from Seq_File_filter.R. This function calls another function Percent_N_apply_Seq.R present in same file. Source the file to load both function into workspace.

 4) The function “Create_consensus_using_UMI.R” identifies sequences with same UMI and a difference of 8 nucleotides or less (PCR repeats), and uses 	them to remove PCR errors and build a consensus based on these sequences.
	Input for this function are :
	i) Files <- list of names of fasta files containing the sequences from each sample( can be one or multiple files).
	ii) name <- Name given to the combined file (all sampled from an organism)
	iii) Forward_primer <- A sequence upstream of the UMI with specific size (For immunoglobulin seq we used the sequence 	"AAGCAGTGGTATCAACGCAGAG"). If a sequence does not have this primer, it is discarded.
	iv) max_mismatch_primer <- maximum mismatch allowed in primer sequence ( For immunoglobulin seq we used 6)
	v) UMI_length <- length of the UMI. It depends on the experimental setup.in this case its 8.
	vi) Mismatch_seq <- mismatch allowed in the sequences with same UMI to be considered as PCR repeats (we used ~2% of sequence length = 8)

 5) All sequences were aligned with the reference sequence using muscle (Multiple Sequence Alignment tool - http://www.drive5.com/muscle/) and edited using the Jalview (multiple sequence alignment editing, visualization and analysis Tool- http://www.jalview.org/). During this process, the first 15 nucleotides of the heavy chain and 3 Nucleotides from the light chain were removed. 

 6) The next step involves replacement of the low quality nucleotides(“N”)with the nucleotides at that position in the reference sequence. The function Remove_N.R does this. The inputs for this function the aligned edited fasta sequences from step 5 and the reference DNA sequence.

 7) The next function Unmutated_remove.R removes all unmutated sequences. The input for this file is the output file from step 6 and the reference DNA sequence.

 8) The function Translate_DNA.R translates the DNA sequences to protein sequence after removing the gaps incorporated during alignment. The input for this file is the output file from step 7. 

 9) The protein sequences are  aligned with the reference sequence using muscle (Multiple Sequence Alignment tool - http://www.drive5.com/muscle/) and edited using the Jalview (multiple sequence alignment editing, visualization and analysis Tool- http://www.jalview.org/). During this step, the non functional protein sequences with stop codon and with frame shift mutations are excluded.

 10) The function “Protein_mutation_stats.R calculated the number of mutations per protein sequence, the number of mutations at each position in all the sequences and the frequency of mutation at each position for all sequences. The input for this function is the aligned filtered protein sequences from the previous step and the reference protein sequence.

 11) DNA_prot_for_R_S.R - This function uses the final filtered functional protein sequences to extract their respective DNA sequences for analysis of replacement and silent mutations. The input files for this function are output DNA sequence fasta file from Unmutated_remove.R and the aligned filtered protein sequences from step 8.

 12)R_S_freq_analysis.R function then calculated the DNA mutation frequency, the frequency of replacement mutation and the frequency of silent mutation for these sequences. The input for this file are the DNA and protein sequences output files from step 10 and the reference protein and DNA sequences.


