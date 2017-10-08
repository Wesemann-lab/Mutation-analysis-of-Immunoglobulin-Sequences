#The function Create_consensus_using_UMI.R finds all sequencins within a file with same UMI, 
#aligns them and create a consensus sequence based on frequency of nucleotides present at each position:
#i) The nucleotide with dominating frequency is kept
#ii) If frequency is equal, the first nucleotide encountered is kept
#iii) If there is an "N" in the position, it is replace by the alternate nucleotide present at the position
#iv) If a single sequence is attached to a UMI, it is left as such


# Input - 
#i) Files <- list of names of fasta files containing the sequences from each sample (can be one or multiple files).
#ii) name <- Name given to the combined file (all sampled from an organism)
#iii) Forward_primer <- A sequence upstream of the UMI with specific size (For immunoglobulin seq we used the sequence "AAGCAGTGGTATCAACGCAGAG"). If a sequence does not have this primer, it is discarded.
#iv) max_mismatch_primer <- maximum mismatch allowed in primer sequence ( For immunoglobulin seq we used 6)
#v) UMI_length <- length of the UMI. It depends on the experimental setup.in this case its 8.
#vi) Mismatch_seq <- mismatch allowed in the sequences with same UMI to be considered as PCR repeats (we used ~2% of sequence length = 8)

# Output
#i)Each file is returned back with a suffix "Unique" with only the unique consensus sequences for each UMI
#ii) A file containing all the uniques sequences combined together
library(Biostrings)
library(seqinr)
library(msa)
Create_consensus_using_UMI <- function(files,name,Forward_primer,max_mismatch_primer,UMI_length,Mismatch_seq){
  #naming the combined file
  file_combined <- paste("Unique_",name, ".fasta", sep = "")
  #loope to run for each file
  for (j in 1 : length(files)){
    #print(j)
    # read the sequnces from the file
    Sequences <- readDNAStringSet(files[j])
    #split the file name to rename each sequence in the file with a number associated with the file
    v_string = strsplit(files[j], "_")
    #count number of sequnces and give number to each sequence along with file name
    count = 1
    #individual output file name - attach "Unique_" to each file name
    file_name <- paste("Unique",files[j], sep = "_")
    # find sequences which have the forward_primer
    m <- vcountPattern(Forward_primer, Sequences,max.mismatch = max_mismatch_primer)
    if(length(which(m==1)) != 0){
      #remove the Forward_primer
      Seq_new <- subseq(Sequences[which(m == 1)], start = length(strsplit(Forward_primer, "")[[1]])+1 )
      #remove duplicate sequnces
      Seq_new <- unique(Seq_new)
      # extract the list of unique UMIs
      Seq_UMI <- unique(subseq(Seq_new, start = 1, end = UMI_length))
      # extract the UMIs of all sequences to compare with the unique UMIs
      Seq_first_eight <- subseq(Seq_new, start = 1, end = UMI_length)
      #for each unique UMI
      for (i in 1: length(Seq_UMI)) {  
        print(i)
        #find all sequences having one UMI
        m <- vcountPattern(Seq_UMI[[i]], Seq_first_eight)
        # extract all the sequences with this UMI and remove the UMI
        mySeq <- subseq((Seq_new[which(m == 1)]), start = UMI_length+1)
        #until all sequences are filtered
        while(length(mySeq) != 0){
          #filter out sequences from mySeq which have a miss match of "Mismatch_seq" or less with the first sequence
          n <- vcountPattern(mySeq[[1]], mySeq ,max.mismatch = Mismatch_seq, fixed = FALSE)
          #if only one sequence with the UMI
          if (length(mySeq) == 1){
            seq = mySeq
            # rename the fasta header with file name and count
            attr(attr(seq,"ranges"),"NAMES") <- paste(name,v_string[[1]][1],Seq_UMI[[i]],count,sep= "_" )
            write.fasta(sequences = seq, names = names(seq), nbchar = 80, file.out = file_name, open = "a")
            write.fasta(sequences = seq, names = names(seq), nbchar = 80, file.out = file_1, open = "a")
            count= count + 1
            rm(seq)
          } else {
            # if no other sequence with mismatch less than "Mismatch_seq"
            if(length(which(n == 1)) < 2){
              seq = mySeq[which(n==1)]
              attr(attr(seq,"ranges"),"NAMES") <- paste(name,v_string[[1]][1],Seq_UMI[[i]],count,sep= "_" )
              write.fasta(sequences = seq, names = names(seq), nbchar = 80, file.out = file_name, open = "a")
              write.fasta(sequences = seq, names = names(seq), nbchar = 80, file.out = file_combined, open = "a")
              count= count + 1
              rm(seq)
            } else {
              # if more than one sequence with mismatch less than "Mismatch_seq" extract these sequences
              Same_UMI <- mySeq[(which(n == 1))]
              # align all sequences 
              myaln <- msa(Same_UMI, "ClustalW")
              #construct consensus matrix
              conMat <- consensusMatrix(myaln)
              conMat <- t(conMat[rowSums(conMat[, -1])>0, ])
              #total = (sum(colSums(conMat)))
              #check frequency and construct the sequence based on the frequency
              consensus = character()
              for (k in 1:dim(conMat)[1]){
                x <- conMat[k,]
                x <- x[order(-x)]
                if(names(x)[1] != "N" && names(x)[1] != "-" )
                {  
                  consensus = rbind(consensus,names(x)[1])
                }  else if (x[2] != 0 && names(x)[2] != "N" && names(x)[2] !=  "-" )
                {
                  consensus = rbind(consensus,names(x)[2])
                } else if (x[3] != 0)
                { 
                  consensus = rbind(consensus,names(x)[3]) 
                } else   {
                  consensus = rbind(consensus, names(x)[1])
                }
                rm(x)   
              }
              #convert the consensus sequence into DNAstring format
              Seq <- DNAStringSet(paste(consensus,collapse=""))
              #write sequence to file
              seq = Same_UMI[1]
              attr(attr(seq,"ranges"),"NAMES") <- paste(name,v_string[[1]][1],Seq_UMI[[i]],count,sep= "_" )
              write.fasta(sequences = Seq, names = names(seq), nbchar = 80, file.out = file_name, open = "a")
              write.fasta(sequences = Seq, names = names(seq), nbchar = 80, file.out = file_combined, open = "a")
              count = count + 1
              rm(seq,Seq,conMat, Same_UMI,k,consensus,myaln)
            }
          }
          # repeat above steps in all sequences until all are processed
          mySeq<-mySeq[(which(n == 0))]
          rm(n)
        }
        rm(m,mySeq)
      }
      rm(Seq_new,Seq_UMI,Seq_first_eight) 
      #print(count)
    }
    rm(Sequences,v_string,count,file_name)
  }
  rm(files,name,Forward_primer,max_mismatch_primer,UMI_length,Mismatch_seq)
} 