library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")



### Run 1 ###

path <- "~/Dropbox/prepro_short/prepro_short/prepro/run1"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(5,5), truncLen=c(228,228),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)
set.seed(101)
errF <- learnErrors(filtFs, multithread=TRUE)
set.seed(101)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

set.seed(101)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
set.seed(101)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
set.seed(101)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
set.seed(101)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))



#save
saveRDS(seqtab, "seqtab_Run1.rds")

#track reads
set.seed(101)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track_r1 <- track
write.table(track_r1, "track_16S_reads_run1.txt")




### Run 2 ###

path <- "~/Dropbox/prepro_short/prepro_short/prepro/run2"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

plotQualityProfile(fnFs[1:12])
plotQualityProfile(fnRs[1:12])



# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(5,5), truncLen=c(228,228),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

set.seed(101)
errF <- learnErrors(filtFs, multithread=TRUE)
set.seed(101)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

set.seed(101)
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
set.seed(101)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

set.seed(101)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
set.seed(101)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

#save
saveRDS(seqtab, "seqtab_Run2.rds")



#track reads
set.seed(101)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track_r2 <- track
write.table(track_r2, "track_16S_reads_run2.txt")


#merge runs and remove chimeras
st1 <- readRDS("seqtab_Run1.rds")
st2 <- readRDS("seqtab_Run2.rds")
set.seed(101)
st.all <- mergeSequenceTables(st1, st2)
set.seed(101)
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)


#assign taxonomy
set.seed(101)
taxa <- assignTaxonomy(seqtab, "silva_138/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
set.seed(101)
taxa <- addSpecies(taxa, "silva_138/silva_species_assignment_v138.1.fa.gz")



#make phyloseq object
mapping_file <- read.table("mapping_rx.txt", header=T)

#match mapping file and otu table and ensure amtch
mapping_file <- mapping_file[match(rownames(seqtab), rownames(mapping_file)),]
all(rownames(mapping_file) == rownames(seqtab))


ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(mapping_file), 
               tax_table(taxa))


