#The function Remove_N.R removes all the low quality nucleotides and replaced them with the residues at that position in the reference sequence

#Input - 
#1) DNA_file -  a fasta file with DNA sequences which have been aligned to the reference sequence using muscle and edited using Jalview
#2)Ref_file - a file containing the reference sequence

#Output -  A fasta file with aligned sequences with no "N" 
library(Biostrings)
library(seqinr)
library(msa)
library(ape)
Remove_N <- function(DNA_file, Ref_file){
  #read the aligned fasta file
  Sequences <- readDNAStringSet(DNA_file)
  #read the reference sequence
  reference <- readDNAStringSet(Ref_file)
  # for each sequence in the DNA_file, check if there is "N in any of the positions and replace it with the nucleotide from the reference sequence at the same position
  for (i in 1:length(Sequences)){
    for(j in 1: length(reference[[1]])){
      if (as.character(Sequences[[i]][j]) == "N"){
        #print(j)
        Sequences[[i]][j] = reference[[1]][j]
      }
    }
  }
  #name of the new file
  name <- paste("No_N",DNA_file, sep = "_")
  #write the Sequences to the new file
  writeXStringSet(Sequences,name)
}

#The function Unmutated_remove.R removes all the unmutated sequences
#Input - 
#1)DNA_file -  a fasta file which is the output from Remove_N.R function
#2)Ref_file - a file containing the reference sequence

#Output -  A fasta file with only unmutated sequences
Unmutated_remove <- function(DNA_file, Ref_file){
  #read the fasta file
  Sequences_DNA <- readDNAStringSet(file)
  #read the reference sequence
  reference <- readDNAStringSet(ref)
  # find sequences which are exactly same as the reference sequence and remove them
  m <- vcountPattern(reference[[1]], Sequences_DNA)
  Sequences_DNA <- Sequences_DNA[which(m==0)]
  #write the rest of the sequnces into a file
  writeXStringSet(Sequences_DNA,file)
}

#The function Translate_DNA.R translates the DNA sequences into protein sequences
#Input - 
#1)DNA_file - The output file from Unmutated_remove.R

#Output -  A fasta file with protein sequences
library(ape)
Translate_DNA<-function(DNA_file)
{
  #Remove gaps to make DNA For translation into protein
  Sequences_DNA <- readDNAStringSet(file)
  Sequences_DNA_prot <- readDNAStringSet(file)
  for (i in 1: length(Sequences_DNA_prot)){
    Sequences_DNA_prot[i] =  gsub("-", "", Sequences_DNA[[i]])
  }
  
  #translate and write to file
  Prot_seq <- trans(as.DNAbin(Sequences_DNA_prot), code = 1, codonstart = 1)
  name <- paste("Protein",file, sep = "_")
  write.dna(Prot_seq, name , format = "fasta", append = "FALSE", nbcol = 80, colsep = "")
}

#The function Protein_mutation_stats.R calculated the mutations in the Protein sequences
#Input - 
#1)Prot_file - The input to this program is the protein sequences aligned to the reference protein sequences (The sequences with stop codons are removed at this step)
#2)Prot_ref - the reference protein sequence
#3)name - sample or filename for the output file

#Output
#1)Prot_mut_per_seq_name <- A text file which gives the number of mutations per sequence
#2)Prot_mut_freq_name <-  A text file which gives the frequency of mutations at each postion for all sequences
#3)Prot_Mutation_matrix_name <- A text file containing a matrix 21 rows and columns equal to the length of the reference sequence. The rows are the amino acids. The entry at each position of the matrix gives us the number of sequences which have mutated at each position and into which amino acid.  

Protein_mutation_stats <- function(Prot_file,Prot_ref,name)
{
  # align and edit the protein and uplaod as edited sequences to calculate the mutation rate for protein
  #read the reference protein sequence file
  ref_prot <- readAAStringSet(Prot_ref)
  #read the aligned and edited protein file
  seq_prot <- readAAStringSet(Prot_file)
  #create consensus Matrix
  conMat_seq <- consensusMatrix(seq_prot)
  colnames(conMat_seq)<-(strsplit(as.character(ref_prot[[1]]), ""))[[1]]
  mutation_mat_prot = conMat_seq
  # if there is no mutation at a specific position, change the count in the consensus matrix to zero
  for ( i in 1: dim(conMat_seq)[2]){
    for( j in  1: dim(conMat_seq)[1]){
      if(colnames(mutation_mat_prot)[i] == rownames(mutation_mat_prot)[j])
      {
        mutation_mat_prot[j,i] = 0
      }
    }
  }
  
  #calculate the number of mutations per sequence
  mutation_per_sequence = matrix(0,nrow =length(seq_prot),ncol = 1 )
  
  for (i in 1:length(seq_prot)){
    mutation_in_seq = 0
    print(i)
    for(j in 1: dim(conMat_seq)[2]){
      if ( seq_prot[[i]][j] != ref_prot[[1]][j]){
        mutation_in_seq = mutation_in_seq + 1
      }
    }
    mutation_per_sequence[i] = mutation_in_seq
  }
  name1 <- paste("Prot_mut_per_seq", name,".txt",sep = "_")
  # write the mutation_per_sequence table to a file
  write.table(mutation_per_sequence, name1, sep = "\t")
  for ( i in 1: dim(conMat_seq)[2])
  {colnames(mutation_mat_prot)[i] = paste(i,colnames(mutation_mat_prot)[i], sep="_" )}
  # calculate the frequency of mutation at each position
  Fre <- (mutation_mat_prot)/length(seq_prot)
  freq_sum <- t(colSums(Fre))
  # write the frequency of mutation at each position to a text file
  name2 <- paste("Prot_mut_freq", name,".txt",sep = "_")
  write.table(freq_sum, name2, sep="\t")
  # write the Mutation_matrix for all sequences to a file
  name3 <- paste("Prot_Mutation_matrix", name,".txt", sep= "_")
  write.table(mutation_mat_prot, name3, sep="\t")
}

#The function DNA_prot_for_R_S.R extracts only functional DNA sequences(No stop codon, no frame shift mutation) based on the aligned and edited protein sequences. These sequences would be used for analysis of silent and replacement mutations
#Input - 
#1)DNA_file - The output DNA sequence fasta file from Unmutated_remove.R
#2)Protein_file - The input prot_file used in the Protein_mutation_stats.R function

#Output -  A new folder is created names R_S, two files are created within it
#1)DNA_R_S.fasta - final filtered DNA sequences for analysis of replcement and silent mutations
#2)Prot_R_S.fasta - final filtered Protein sequences for analysis of replcement and silent mutations

#note: this function requires for the names of the DNA sequences to match with their respective protein sequences

#extract the dna sequence of functional proteins only and write them for r/s analysis    
DNA_prot_for_R_S <- function(DNA_file, Protein_file){
  if(!file.exists("R_S"))
  { dir.create("R_S")}
  
  #read the fasta file for DNA sequences
  Seq_DNA <- readDNAStringSet(DNA_file)
  #read the fasta file for Protein sequences
  Seq_Prot <- readAAStringSet(Protein_file)
  Seq_names<- names(Seq_DNA)
  # if the names of two sequences in the two data sets match, they are written to respective DNA and protein fasta files
  for(j in 1 : length(Seq_Prot)) {
    m <- vcountPattern(names(Seq_Prot[j]), Seq_names)
    write.fasta(sequences = Seq_Prot[j], names = names(Seq_Prot[j]), nbchar = 80, file.out = "R_S/Prot_R_S.fasta", open = "a")
    write.fasta(sequences = Seq_DNA[which(m == 1)], names = names(Seq_DNA[which(m == 1)]), nbchar = 80, file.out = "R_S/DNA_R_S.fasta", open = "a")
  }
}

#The function R_S_freq_analysis.R calculated the frequencies of total mutation in the DNA sequences, the frequency of replacement mutation and the frequency of silent mutation
#Input - 
#1)DNA_file - The output DNA sequence fasta file from DNA_prot_for_R_S.R function
#2)Protein_file - The ioutput Protein sequence fasta file from DNA_prot_for_R_S.R function
#3)DNA_ref - reference DNA sequence
#4)Prot_ref - reference Protein sequence

#Output -  
#1)R_S_dna.txt - The first row of this table gives us the total number of mutations at each position in all DNA sequences, The second row gives the number of sequences in which the mutation caused change in amino-acid at the position while the third row gives the number of sequences with silent mutations at respective positions.  
#2)R_S_freq.txt - The first row of this table gives us the frequency of mutations at each position in all DNA sequences, The second row gives the frequency of mutation which cause change in amino-acid at the position while the third row gives the frequency of sequences with silent mutations at respective positions.

R_S_freq_analysis <- function(DNA_file, Prot_file, DNA_ref, Prot_ref)     {
  Prot_ref <- readAAStringSet(Prot_ref)
  Prot_seq <- readAAStringSet(Prot_file)
  #calculate total number of sequences
  total_seq <- length(Prot_seq)
  DNA_ref <- readDNAStringSet(DNA_ref)
  DNA_seq <- readDNAStringSet(DNA_file)
  
  Prot_refer <- (strsplit(as.character(Prot_ref[[1]]), ""))[[1]]
  ref_dna_vec <- (strsplit(as.character(DNA_ref[[1]]), ""))[[1]]
  
  #create a matrix twith three rows and number of columns equal to the length of DNA sequence. This matrix will store the Total mutation, Replacement and Silent mutation number.
  R_S <- matrix(0,nrow = 3,ncol = length(DNA_ref[[1]]))
  colnames(R_S) <- (strsplit(as.character(DNA_ref[[1]]), "" ))[[1]]
  rownames(R_S) <- c("DNA_mutation","Replacement","Silent")
  #for each sequence, check if there was a mutation, if yes, check if the mutation lead to replacement of amino acid or not
  for (i in 1: total_seq) 
  { #print ("This is i")
    #print(i)
    DNA <- strsplit(as.character(DNA_seq[[i]]), "")[[1]]
    Protein <-   strsplit(as.character(Prot_seq[[i]]), "")[[1]]
    for(k in 1:length(DNA_ref[[1]]) ){
      if(DNA[k] != ref_dna_vec[k]){
        R_S["DNA_mutation",k] =  R_S["DNA_mutation",k] + 1
        R_S["Silent",k] =   R_S["Silent",k] + 1
        if(k%%3 == 0)
        {l = k%/%3} else {l = (k%/%3)+1}
        if(Protein[l] != Prot_refer[l]) {
          R_S["Replacement",k] =   R_S["Replacement",k] + 1
          R_S["Silent",k] = R_S["Silent",k] - 1
        }
        #if(Protein[l] == "-"){  R_S["Silent",k] = "D"}
      }
    }
  }
  #write the table to a file
  write.table(R_S,"R_S_dna.txt", sep = "\t")
  #calculate the frequency for the above table 
  Frequency_R_S <- matrix(0,nrow = 3,ncol = length(DNA_ref[[1]]))
  colnames(Frequency_R_S) <- (strsplit(as.character(DNA_ref[[1]]), "" ))[[1]]
  rownames(Frequency_R_S) <- c("DNA_mutation_Freq","Replacement_Freq","Silent_freq")
  
  for( i in 1 : dim(Frequency_R_S)[2])
  {
    Frequency_R_S[1,i] = R_S[1,i]/total_seq
    Frequency_R_S[2,i] =R_S[2,i]/total_seq
    Frequency_R_S[3,i] = R_S[3,i]/total_seq
    #Frequency_R_S[4,i] = length(which(Sil_rec[,i] == "D"))/total_seq
  }
  #write the frequency file
  write.table(Frequency_R_S,"R_S_freq.txt", sep = "\t")
}
