#!/usr/bin/env Rscript

# post-processing salmon run
# post-process-salmon.R
# this command is called from a bash script running in a particular location
# relative to the main directory of data
#setwd("/Users/mjk/Desktop/Tresorit_iOS/projects/RNA-Seq/upstreamRNASeq/workflow")

# I need tc -> gene mapping file which will be in stringtie_merged.gtf file

# I need TMM function

# I need flags for doing gene-level quant and for using TMM function or not


require(data.table)
require(magrittr)
require(dplyr)
require(plyr)


###############
# OPTIONAL ARGS:
args = commandArgs(trailingOnly=TRUE)
geneLevel = args[1]

 geneLevel = TRUE

###############
#  FUNCIONS:

initializeDF <- function(df, colName1=NA){
   j = ncol(df) - 1
   colClasses <- c("character", rep("numeric", j))
   if(is.na(colName1)){col.names <- colnames(df)}
   else{col.names <- c(colName1, colnames(df)[2:ncol(df)])}
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
   # if not, then maps it to the gene id.
   lookup <- data.frame(tc="", gene="")
   for(line in annotID){
      flds <- unlist(strsplit(line, ";"))
      gene_id <- gsub("\"$","",gsub("gene_id \"","",flds[1]))
      tc_id <- gsub("\"$","",gsub(" transcript_id \"","",flds[2]))
      if(!is.na(flds[3])){
         gene_name <- gsub("\"$","",gsub(" gene_name \"","",flds[3]))
         lookup <- rbind.fill(lookup, data.frame(tc=tc_id, gene=gene_name))
         #print(c(gene_id, tc_id, gene_name))
      }
      else{
         #print(c(gene_id, tc_id))
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
df <- read.csv(quantFiles[1], header=TRUE, stringsAsFactors = FALSE, sep="\t")[c("Name", "TPM")]
colnames(df) <- c("Name", grabName(quantFiles[1]))

# extend df
for(file in quantFiles[2:length(quantFiles)]){
   quant <- read.csv(file, header=TRUE, stringsAsFactors = FALSE, sep="\t")[c("Name", "TPM")]
   if(any(quant$Name != df$Name)) stop("Transcript Names don't align")
   else{
      df <- cbind(df, quant$TPM)
      colnames(df)[dim(df)[2]] <- grabName(file)
   }
}



# gene-level?
# mapping transcript_id to (a) gene_name or (b) gene_id if (a) doesn't exist
if(!is.na(geneLevel)){
   mergedGTF <- fread("stringtie_merged.gtf", sep="\t", skip=2)
   mergedGTF <- mergedGTF %>% filter(V3=="transcript")
   tc2gene <- processFields(mergedGTF$V9)
   GeneData <- initializeDF(df, colName1 = "Gene")
   #GeneData <- data.frame(Gene=character(), ERR188044_chrX=numeric(), 
   #                       ERR188104_chrX=numeric())
   # for each gene, grab all TPMs for various Tcs and sum, creating DF
   for(Gene in unique(tc2gene$gene)){
      if(Gene==""){next}
      print(Gene)
      rows <- tc2gene %>% filter(gene==Gene)
      Tc_data <- df %>% filter(Name %in% rows$tc)
      print(Tc_data)
      print(c("dim Tc_data", dim(Tc_data)))
      
      #print(rows)
      #Tc_data <- initializeDF(df)
      #<- data.frame(Name=NA, ERR188044_chrX=NA, ERR188104_chrX=NA)
### AUTOMATE PREVIOUS STATEMENT OF DEFINING COLNAMES
      # for(tc in rows$tc){
      #    # grab data from 'df' that match tc
      #    tmp <- df %>% filter(Name==tc)
      #    names(tmp) <- names(df)
      #    Tc_data <- rbind(Tc_data, tmp)
      #    print(head(Tc_data))
      #    print(class(Tc_data[1,2]))
      #    print(class(Tc_data[1,1]))
      #    #Tc_data <- rbind(Tc_data, df %>% filter(Name==tc))
      # }
      # skip lines that are empty:
      #if(nrow(Tc_data) > 1){
      gd <- c(Gene, apply(Tc_data[,2:ncol(Tc_data)], 2, sum))
      names(gd)[1] <- "Gene"
      print(c("dimensions gd", dim(gd)))
      #colnames(gd) <- colnames(GeneData)
      print(head(gd))
      #print(class(gd[1,2]))
      GeneData <- rbind(GeneData, gd, stringsAsFactors=F)
     # }
      #else{
         # print(c(data, "***"))
      #}
   }
}
GeneDataDF <- cbind(GeneData[,1], data.frame(apply(GeneData[,2:dim(GeneData)[2]],
                                                   2, as.numeric)))
# filter out NA's
GeneDataDF <- GeneDataDF %>% filter(!is.na(GeneDataDF[,2]))


if(is.na(geneLevel)){                       # return Tc-level df
   
}


# test code
# I want to automatically create empty dataframes with defined column names
# the first name will be user-inputted and all others will be read, ultimately
# from salmonQfiles.txt.  The first class will be character and all others will
# be numeric
colnames(df)
#catString = paste0(colnames(df)[1], "=character()")
#for(name in colnames(df)[2:dim(df)[2]]){
#   catString <- paste0(catString, ", ", name, "=numeric()")
#}
#catString <- paste0(catString, ", stringsAsFactors=FALSE")
#df2 <- data.frame(eval(catString))




tmp <- list("Sam", 18, 14)
names(tmp) <- names(df5)
rbind(df5, tmp)



