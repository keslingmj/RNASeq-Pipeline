#!/usr/bin/env Rscript

# post-process-salmon.R
# this command is called from a bash script running in a particular location
# relative to the main directory of data

# Michael Kesling, 2019-11-27


require(data.table)
require(magrittr)
require(dplyr)
require(plyr)
library(tibble)


###############
# OPTIONAL ARGS:
args = commandArgs(trailingOnly=TRUE)
geneLevel = args[1]

 geneLevel = TRUE

###############
#  FUNCIONS:

initializeDF <- function(df){  
   j = ncol(df) - 1
   colClasses <- rep("numeric", ncol(df))
   col.names <- colnames(df)
   tmpDF <- read.table(text="",
                       colClasses=colClasses,
                       col.names = col.names,
                       stringsAsFactors = FALSE)
   return(tmpDF)
}


grabName <- function(string){
   return(gsub(".quant/quant.sf","", gsub("salmonQuants/", "", string)))
}

processFields <- function(annotID){
   # maps each transcript ID to a gene name, if available.
   # if not, then maps it to the gene id. Uses stringtie_merged.gtf column 9 as input
   lookup <- data.frame(tc="", gene="")
   for(line in annotID){
      flds <- unlist(strsplit(line, ";"))
      gene_id <- gsub("\"$","",gsub("gene_id \"","",flds[1]))
      tc_id <- gsub("\"$","",gsub(" transcript_id \"","",flds[2]))
      if(!is.na(flds[3])){
         gene_name <- gsub("\"$","",gsub(" gene_name \"","",flds[3]))
         lookup <- rbind.fill(lookup, data.frame(tc=tc_id, gene=gene_name))
      }
      else{
         lookup <- rbind.fill(lookup, data.frame(tc=tc_id, gene=gene_id))
      }
   }
   return(lookup)
}


##########
# MAIN
##########

quantFiles <- readLines("salmonQfiles.txt")

# initialize df with first quantFile:
df <- read.table(quantFiles[1], header=TRUE, stringsAsFactors = FALSE, sep="\t", row.names=1)["TPM"]
colnames(df) <- grabName(quantFiles[1])

# extend df
for(file in quantFiles[2:length(quantFiles)]){
   quant <- read.csv(file, header=TRUE, stringsAsFactors = FALSE, sep="\t", row.names=1)["TPM"]
   if(any(rownames(quant) != rownames(df))) stop("Transcript Names don't align")
   else{
      df <- cbind(df, quant$TPM)
      colnames(df)[ncol(df)] <- grabName(file)
   }
}

# gene-level summation of TPMs here
# mapping transcript_id to (a) gene_name or (b) gene_id if (a) doesn't exist
if(!is.na(geneLevel)){
   mergedGTF <- fread("stringtie_merged.gtf", sep="\t", skip=2)
   mergedGTF <- mergedGTF %>% filter(V3=="transcript")
   tc2gene <- processFields(mergedGTF$V9)
   GeneData <- initializeDF(df)
   
   # for each gene, grab all TPMs for various Tcs and sum, creating gd DF
   # and add to growing GeneData DF
   geneNames <- c()
   for(Gene in unique(tc2gene$gene)){
      if(Gene==""){next}
      rows <- tc2gene %>% filter(gene==Gene)
      Tc_data <- df %>% rownames_to_column('Name') %>%
         filter(Name %in% rows$tc) %>% column_to_rownames('Name')
      
      gd <- apply(Tc_data, 2, sum)
      names(gd) <- colnames(Tc_data)
      
      GeneData <- rbindlist(list(GeneData, as.list(gd)))          # (data.table)
      geneNames <- c(geneNames, Gene)
   }
   # convert GeneData to dataframe, then add rownames
   GeneData <- as.data.frame(GeneData)
   rownames(GeneData) <- geneNames
}

###########
# write file to file system
write.table(GeneData, "salmonTPMmerged.txt", sep="\t", row.names = T,
            col.names = NA)         # NA needed for offsetting col names in unix



if(is.na(geneLevel)){                       # return Tc-level df
   # code not yet written
   return("MUST USE 'TRUE' option, as only gene-level expression currently supported")
}