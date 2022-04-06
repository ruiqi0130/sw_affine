#!/usr/bin/env Rscript

### Usage: Rscript --vanilla hw1.R <input file> <score file>
### Example: Rscript --vanilla hw1.R input.txt blosum62.txt
### Note: Smith-Waterman Algorithm

### This is one way to read in arguments in R
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
stop("At least two arguments must be supplied (inputFile, scoreFile).n", call.=FALSE) } else if (length(args)>=2) {
    # default gap penalties
    args[3] = -2
    args[4] = -1 }

## Specifying author and email
p <- c(person("Ruiqi", "Li", role = "aut", email = "ruiqi.li@yale.edu"))

## Implement your Smith-Waterman Algorithm
runSW <- function(inputFile = "./input.txt", scoreFile = "./blosum62.txt", openGap = -2, extGap = -1) {
    ### Read input
    str1 = read.table(inputFile, header = FALSE)[1,]
    str2 = read.table(inputFile,header = FALSE)[2,]
    
    str1 = as.character(unlist(strsplit(str1,split="")))
    str2 = as.character(unlist(strsplit(str2,split="")))
    
    score_matrix = read.table(scoreFile, header = TRUE)
    
    ### Calculate similarity matrix
    sim_matrix = matrix(0, ncol = length(str1)+1, nrow = length(str2)+1)
    colnames(sim_matrix) = c("\t",str1)
    rownames(sim_matrix) = c(" ",str2)
    
    for(i in 2:nrow(sim_matrix)){
        for(j in 2:ncol(sim_matrix)){
            M_diag = sim_matrix[i-1,j-1] + score_matrix[rownames(sim_matrix)[i],
                                                        colnames(sim_matrix)[j]]
            M_up = max( unlist(sapply(1:(i-1),
                                       function(k){sim_matrix[i-k,j] + 
                                               openGap + (k - 1) * extGap})) )
            M_left = max( unlist(sapply(1:(j-1),
                                      function(k){sim_matrix[i,j-k] + 
                                              openGap + (k - 1) * extGap})) )
            sim_matrix[i,j] = max(M_diag,M_up,M_left,0)
        }
    }
    
    ### Trace Back
    align_score = max(sim_matrix)
    row_start = which(sim_matrix == align_score, arr.ind = TRUE)[1]
    col_start = which(sim_matrix == align_score, arr.ind = TRUE)[2]
    
    match_str1 = NULL
    match_str2 = NULL 
    match_align = NULL
    r = row_start; c = col_start

    while(sim_matrix[r,c]!=0){
        # update match_align

        if(sim_matrix[r-1,c-1] + score_matrix[rownames(sim_matrix)[r],
                                              colnames(sim_matrix)[c]] == sim_matrix[r,c]){
            match_str1 = c(colnames(sim_matrix)[c],match_str1)
            match_str2 = c(rownames(sim_matrix)[r],match_str2)
            if(colnames(sim_matrix)[c] == rownames(sim_matrix)[r]){
                match_align = c("|",match_align)
            }else{
                match_align = c(" ",match_align)
            }
            r = r-1; c = c-1
        }
        else if(
            any(sapply(1:(r-1),function(k){
                sim_matrix[r-k,c] +
                    openGap + (k - 1) * extGap == sim_matrix[r,c]}))
        ) {
            k = which(sapply(1:(r-1),function(k){
                sim_matrix[r-k,c] + 
                    openGap + (k - 1) * extGap == sim_matrix[r,c]}) )
            match_str1 = c(rep("-",k),match_str1)
            match_str2 = c(rownames(sim_matrix)[(r-k+1):r],match_str2)
            match_align = c(rep(" ",k),match_align)
            r = r-k; c = c
        }
        # else if(
        #     any(sapply(1:(c-1),function(k){
        #         sim_matrix[r,c-k] + 
        #             openGap + (k - 1) * extGap == sim_matrix[r,c]})) 
        # ) 
        else{
            k = which(sapply(1:(c-1),function(k){
                sim_matrix[r,c-k] + 
                    openGap + (k - 1) * extGap == sim_matrix[r,c]}) )
            match_str1 = c(colnames(sim_matrix)[(c-k+1):c],match_str1)
            match_str2 = c(rep("-",k),match_str2)
            match_align = c(rep(" ",k),match_align)
            r = r; c = c-k
        }
    }
         
    
    #### joint the rest string
    ##### add the rest part
    substr1_l = str1[c(1:c-1)]
    substr2_l = str2[c(1:r-1)]
    if(col_start <= length(str1)){
        substr1_r = str1[c(col_start:length(str1))]
    }else{substr1_r = NULL}
    if(row_start <= length(str2)){
        substr2_r = str2[c(row_start:length(str2))]
    }else{substr2_r = NULL}
    
    ##### add blank
    if(length(substr1_l) > length(substr2_l)){
        substr2_l = c(rep(" ",length(substr1_l)),substr2_l)
    }
    substr1_l = c( rep(" ",max(length(substr2_l)-length(substr1_l),0)) , substr1_l)
    substr2_l = c( rep(" ",max(length(substr1_l)-length(substr2_l),0)) , substr2_l)
    substr1_r = c( rep(" ",max(length(substr2_r)-length(substr1_r),0)) , substr1_r)
    substr2_r = c( rep(" ",max(length(substr1_r)-length(substr2_r),0)) , substr2_r)
    
    match_str1 = c(substr1_l,"(",match_str1,")",substr1_r)
    match_str2 = c(substr2_l,"(",match_str2,")",substr2_r)
    match_align = c( rep(" ",length(substr1_l)+1),match_align,
                     rep(" ",length(substr1_r)+1) )

    ### Write to output file
    str_input = c(paste(str1,collapse = ""),
                  paste(str2,collapse = ""))
    str_output = c(paste(match_str1,collapse = ""),
                   paste(match_align,collapse = ""),
                   paste(match_str2,collapse = ""))
    outputFile = "test.txt"
    cat("-----------\n|Sequences|\n-----------\n", file=outputFile)
    cat("sequence1\n",str_input[1],"\nsequence2\n",str_input[2],"\n", sep = "",
        file = outputFile, append = TRUE)
    cat("--------------\n|Score Matrix|\n--------------\n", 
        file = outputFile, append = TRUE)
    write.table(sim_matrix, sep = "\t", quote = FALSE, col.names=TRUE,
                file = outputFile, append = TRUE)
    cat("----------------------\n|Best Local Alignment|\n----------------------\n",
        file = outputFile, append = TRUE)
    cat("Alignment Score:",align_score,"\nAlignment Results:\n",
        str_output[1],"\n",str_output[2],"\n",str_output[3],"\n",sep="",
        file = outputFile, append = TRUE)
}


## Run the main function and generate results
runSW(inputFile=args[1], scoreFile=args[2], openGap=args[3], extGap=args[4])