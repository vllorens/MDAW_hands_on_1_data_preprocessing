# Microbiome data analysis workshop
# Hands on #1: Raw data to OTU table
# Instructor: Veronica Llorens-Rico, PhD
# Date: 20/04/2021

#### PART 1: Introduction ####
# Short pptx presentation about the following topics:
    
# - High-throughput 16S amplicon sequencing
# - First preprocessing steps from the raw data: removal of primers, demultiplexing
# - Differences between OTUs and ASVs
# - Examples of pipelines to preprocess the data


#### PART 2: from the *.fastq files to the OTU table: using the DADA2 pipeline ####
# This part largely follows the DADA2 tutorial from https://benjjneb.github.io/dada2/tutorial.html
# But with some modifications adapted to our own testing dataset

# The data that we will use was originally generated for the following publication:
# Quantitative microbiome profiling links gut community variation to microbial load 
# (https://www.nature.com/articles/nature24460)
# In this article, the gut microbiome of 40 healthy adults was profiled via 16S sequencing, 
# and absolute microbial loads per gram of stool were determined via flow cytometry. 
# Here, we will reprocess the data using DADA2 for 20 of those healthy adults. 
# Later in the course, I will also give a seminar on data normalization and 
# why determining microbial loads is important in microbiome research.

# To download the data, follow these steps:
    
# 1. Go to https://www.ebi.ac.uk/ena/browser/view/PRJEB21504
# 2. From the "Submitted FTP" column, select samples SC01 to SC20 (both R1 and R2 fastq files), you will need to go through the first 4 pages
# 3. Check that the option "Download files as zip" is active
# 4. Click on "Download selected files"
# 5. Move the "ena_files.zip" to the `data/` folder from the course and unzip.

## 2.1. Package installation and checks
packages_cran <- c("BiocManager", "ggplot2", "ggpubr")
for(pkg in packages_cran){
    if(!pkg %in% installed.packages())
        install.packages(pkg)
}

packages_bioconductor <- c("dada2", "phyloseq")
for(pkg in packages_bioconductor){
    if(!pkg %in% installed.packages())
        BiocManager::install(pkg, update = F)
}

# make sure DADA2 package is at least v1.12 - otherwise update the package
if(packageVersion("dada2")<"1.12"){
    install.packages("dada2")
}


## 2.2. Load required packages
library(dada2)

## 2.3. Check that the you can find the files in the data directory

path <- "data/" # When working with your own data, change the path to the directory containing your demultiplexed R1 and R2 fastq files
fastqseqs <- list.files(path, pattern = ".fastq.gz") # take only fastq files, ignore the zip file
fastqseqs

# You should see a list of 40 files in the output

fastqseqs <- sort(fastqseqs) # Sort ensures forward/reverse reads are in same order (they should already be)

## 2.4. Separate forward (R1) from reverse (R2) files and select sample names

# Make sure that R1 is for forward read and R2 for reverse! 
# In some sequencing protocols, the reads are reversed.
fnFs <- fastqseqs[grepl(".R1.fastq.gz", fastqseqs)] # Just select forward read files
fnRs <- fastqseqs[grepl(".R2.fastq.gz", fastqseqs)] # Just select reverse read files

# Get sample names from the first part of the forward read filenames
sample.names <- sapply(strsplit(fnFs, ".R1.fastq.gz"), `[`, 1) ## check if is 1 or 2!

# For the processing, get the full path to the files, not just the file names
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

## 2.5. Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2]) 

# The quality of the reverse read decays after ~130 bases, so we will trim it. 
# This is normal for Illumina, so there's nothing to worry about

## 2.6. Filter and trim sequences

# Place filtered files in filtered/ subdirectory in the data folder
filt_path <- file.path(path, "filtered") 

# Create names for the filtered files
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter - this takes ~2 min in a Macbook from 2015
out <- filterAndTrim(fnFs, # original forward files
                     filtFs, # filtered forward files
                     fnRs, # original reverse files
                     filtRs, # filtered reverse files
                     truncLen=c(200,140), # these settings are good for MiSeq (for our primer constructs)
                     trimLeft=c(30, 30), # trim the first bases (if low quality or if primers have not been removed)
                     maxN=0, # discard sequences with >0 Ns after truncation
                     maxEE=c(2,2), # discard reads with more than 2 Expected Errors (calculated from the quality score)
                     truncQ=11, # truncate reads when quality decreases below 11
                     rm.phix=TRUE, # removes traces of PhiX phage, used as internal control for sequencing
                     compress=TRUE, # compress the output
                     multithread=TRUE) # change to F if your computer does not have several cores you can use
head(out) 

# After filtering, we can check the quality of the filtered files to verify that we have solved the issues
plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])

## 2.7. Learn error rates
# NOTE: for speed reasons, the argument nbases for the learnErrors functions is quite low in this tutorial, 
# for real cases, it's better use the default 1e8 
# With nbases = 1e7; this takes ~2min to run in a Macbook from 2015
set.seed(12345)
# Learn forward error rates 
errF <- learnErrors(filtFs, nbases=1e7, multithread=TRUE) 
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e7,multithread=TRUE) 

# plot error rates
plotErrors(errF, nominalQ=TRUE)


## 2.9. Sample inference 
# In this tutorial, we skip the step of dereplication, 
# because as of DADA2 v1.12, it is performed within the main dada function
# This part takes ~8 minutes in a macbook from 2015
dadaFs <- dada(filtFs, # the function will use the filtered files
               err=errF, # we pass on the calculated error rates
               multithread=TRUE # change if your computer does not allow multithreading
               )
dadaRs <- dada(filtRs, # the function will use the filtered files
               err=errR, # we pass on the calculated error rates
               multithread=TRUE # change if your computer does not allow multithreading
               )


## 2.10. Merge forward and reverse reads after dereplication and denoising
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
names(mergers) <- sample.names
head(mergers[[1]])


## 2.11. Construct sequence (ASV) table and remove chimeras
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

## 2.12. Last checks: track reads throughout the DADA2 pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


#### PART 3. Assign taxonomy to the ASV table ####
# Different databases of taxonomic assignments of 16S sequences exist, here, we will use the RDP database
# Assign taxonomy (~5 min on a macbook from 2015)
date()
taxHS <- assignTaxonomy(seqtab.nochim, # the ASV table
                        refFasta = "data/rdp_train_set_16.fa.gz",
                        tryRC = T,
                        multithread=TRUE)
date()

# Add species
taxHS <- addSpecies(taxHS, "data/rdp_species_assignment_16.fa.gz", tryRC = T)
date()

colnames(taxHS) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
unname(head(taxHS))
unname(tail(taxHS))

# Write to disk
write.table(track, file = "output/track_reads.tsv", quote=FALSE)
write.table(seqtab.nochim, file = "output/sequence_table_SV.tsv", quote=FALSE)
write.table(taxHS, file = "output/taxonomy_SV.tsv", quote=FALSE)


#### PART 4. Transfer data to phyloseq ####

library(phyloseq)

# Phyloseq objects contain all the information about the 16S data, structured in 3 parts:
# - otu_table (the table containing the frequencies of all OTUs or, in our case, ASVs)
# - tax_table (the table with the taxonomic annotation)
# - sample_data (metadata of the samples)

# Here, we don't have metadata, but we have the other two components, 
# and we can build a phyloseq object with them
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F),
               tax_table(taxHS))

ps

# Phyloseq allows for the manipulation and further analyses of microbiome data.
# As an example, here we will plot the top 15 taxa of our sample (genus level)
# Further info on the phyloseq package will come on further tutorials

# first, we agglomerate the data at the genus level 
ps_genus <- tax_glom(ps, "Genus")

# select the top 15 genera 
top15 <- names(sort(taxa_sums(ps_genus), decreasing=TRUE))[1:15]

# calculate relative abundances
ps_genus <- transform_sample_counts(ps_genus, function(ASV) ASV/sum(ASV))

# group everything that's not in the top 10 as NA and plot
tax_table(ps_genus)[!rownames(tax_table(ps_genus)) %in% top10,"Genus"] <- NA
plot_bar(ps_genus, fill="Genus", title = "Top 15 taxa") 







