library(plotrix)
library(seqinr)
library(dplyr)
library(tidyr)


### QUESTION 1: GENOME ASSEMBLY
statsnocutoff <- read.delim("statsnocutoff.txt")
hist(statsnocutoff$short1_cov, breaks=10000, xlim=c(0,50))
weighted.hist(statsnocutoff$short1_cov, statsnocutoff$lgth, breaks=0:50, ylab="Frequency")

statscov8exp21 <- read.delim("statscov8exp21.txt")
hist(statscov8exp21$short1_cov, breaks=10000, xlim=c(0,50))
weighted.hist(statscov8exp21$short1_cov, statscov8exp21$lgth, breaks=0:50)

statscov18exp21 <- read.delim("statscov18exp21.txt")
hist(statscov18exp21$short1_cov, breaks=10000, xlim=c(0,50))
weighted.hist(statscov18exp21$short1_cov, statscov18exp21$lgth, breaks=0:50)

statsauto <- read.delim("statsauto.txt")
hist(statsauto$short1_cov, breaks=10000, xlim=c(0,50))
weighted.hist(statsauto$short1_cov, statsauto$lgth, breaks=0:50)

  
#N50 ~400,000 no.nodes ~300? used ~150,000

#high coverage:
  #repeats
  #PCR artefacts 


contigscov18exp21 <- read.fasta("contigscov18exp21.fa")

find.length <- function(string){
  string <- paste0(string, collapse="")
  length <- nchar(string)
  #result <- c(length,string)
  #result
  length - 30 #fixed artefact number
}

lengths.contigs <- sapply(contigscov18exp21, FUN=find.length )
sorted.lengths.contigs <- order(lengths.contigs, decreasing = T)
sorted.lengths.contigs[1:10]

max.node <- lengths.contigs[sorted.lengths.contigs[1]]
max.node

top.ten.nodes <- lengths.contigs[sorted.lengths.contigs[1:10]]
top.ten.nodes

write.fasta(paste0(contigscov18exp21$NODE_69_length_6849_cov_24.023069 , 
                   collapse=""), names = "NODE_69_length_6849_cov_24.023069", 
            file.out="longcontig2.fa")



### QUESTION 2: SPECIES IDENTITY AND GENE ALIGNMENT W EXONERATE
# exonerate Buchnera_aphidicola_bcc.ASM9096v1.cdna.all.fa Escherichia_coli_str_k_12_substr_mds42.GCA_000350185.1.23.cdna.all.names.fa --score 100 --ryo ">> %qi \t %ti \t %pi \t %s \t %ql \n" --showalignment no --showvulgar no  > exonerate_score100
# grep ">>" exonerate_output | awk '{if((NR)>1) print $2,$3,$4,$4,$6}' > test.txt 


exonerate.output <- read.delim("exonerate_output.txt", header = F, sep="\t")
foo <- function(input){
  split.string <- strsplit(as.character(input), split=" ")
  split.string
}
exonerate.output.list <- lapply(exonerate.output,FUN = foo)
exonerate.output.df <- data.frame(matrix(unlist(exonerate.output.list), nrow=520, byrow=T),stringsAsFactors=FALSE)
colnames(exonerate.output.df) <- c("qi", "ti", "percentage.similarity", "raw.score", "query.length")

exonerate.output.df <- exonerate.output.df[complete.cases(exonerate.output.df), ]
exonerate.output.df <- exonerate.output.df[-grepl("BCc", exonerate.output.df$qi), ]
dim(exonerate.output.df)

exonerate.output.df$raw.score <- as.numeric(exonerate.output.df$raw.score)

hist(exonerate.output.df$raw.score, breaks=1000,xlim=c(0,200), xlab="Raw score", main="Histogram of raw score")

#Histogram shows that a threshold of 125 needs to be set for the raw score


exonerate.output.thresholded <- exonerate.output.df[c(as.numeric(exonerate.output.df$raw.score) > 125), ]
rownames(exonerate.output.thresholded) <- NULL
exonerate.output.thresholded


### QUESTION 3: DIFFERENCES IN GENE CONTENT BETWEEN SPECIES

Buchnera <- read.fasta("Buchnera_aphidicola_bcc.ASM9096v1.cdna.all.fa")
E.Coli <- read.fasta("Escherichia_coli_str_k_12_substr_mds42.GCA_000350185.1.23.cdna.all.names.fa")

#lengths of genomes 
length(Buchnera) #365
length(E.Coli) #3635


#number of matched genes
dim(exonerate.output.thresholded) #117 matched genes with raw score >125

#number of unique Buchenara genes mapped onto E.Coli
unique.buchenara <- unique(exonerate.output.thresholded$qi)
length(unique.buchenara)  #86 unique Buchnara genes

#percentage of Buchenara in E.Coli:
(length(unique.buchenara))/length(Buchnera)*100 #23.5%

#percentage of E.Coli in Buchenara:
unique.ecoli <- unique(exonerate.output.thresholded$ti)
length(unique.ecoli) /length(E.Coli) *100 #2.5%


# Duplicated Buchnara genes
duplicated.qi <- duplicated(exonerate.output.thresholded$qi)
sum(duplicated.qi) #31 duplicates genes
#Unique Buchnara genes that are acting as many-to-one matching
length(unique(exonerate.output.thresholded$qi[duplicated.qi])) #21 unique genes


#Duplicated E.Coli genes
duplicated.ti <- duplicated(exonerate.output.thresholded$ti)
sum(duplicated.ti) #25 duplicating genes
length(unique(exonerate.output.thresholded$ti[duplicated.ti])) #19 unique genes




#Histogram of percentage id
hist(exonerate.score100$`%pi`, breaks=100, xlab="Percentage Similarity", main="Histogram of percentage similarity of matched genes", ylim=c(0,50), xlim=c(0,100))


unique.ecoli.separated <- paste0(unique.ecoli, collapse="", sep='\n')
cat(unique.ecoli.separated)
length(unique.ecoli)
unique.ecoli


### Gene ontology analysis

go.results <- read.csv("GOanalysis.txt")
print(go.results[9,1])
    