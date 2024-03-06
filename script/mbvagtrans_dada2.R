# Script for 16S Analysis of MBVag Trans project (2023)
# Raw data are 16S sequencing results (Fastq files) in project folder 22010/22010_RawData/, this folder is not included in git
library(dada2)
library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(phangorn)

# load fastq files 
path = "22010/22010_RawData/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# run dada2
library(dada2); packageVersion("dada2")
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path("intermediate", "fq_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("intermediate", "fq_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(17,21),
                     truncLen=c(245,245),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrF = plotErrors(errF, nominalQ=TRUE)
plotErrR = plotErrors(errR, nominalQ=TRUE)

saveRDS(errF, file = "intermediate/dada2/errF.rds")
saveRDS(errR, file = "intermediate/dada2/errR.rds")

ggsave(plotErrF, filename = "results/dada2/errF.png")
ggsave(plotErrR, filename = "results/dada2/errR.png")

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

saveRDS(dadaFs, file = "intermediate/dada2/dadaFs.rds")
saveRDS(dadaRs, file = "intermediate/dada2/dadaRs.rds")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
saveRDS(seqtab, file = "results/dada2/seqtab.rds")

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # overall frequence of non-chimeric reads is 87.5 %
saveRDS(seqtab.nochim, file = "results/dada2/seqtab_nochim.rds")

# track reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track_long = as.data.frame(track) %>%
  rownames_to_column("sample") %>%
  tidyr::pivot_longer(input:nonchim, names_to = "step", values_to = "number") %>%
  mutate(step = factor(step, levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")))

trackplot = ggplot(track_long, aes(x = step, y = number)) +
  geom_boxplot() +
  geom_jitter()+
  theme_bw(base_size = 14) +
  scale_y_log10()

write_excel_csv2(as.data.frame(track), file = "results/dada2/track.csv")
ggsave(trackplot, filename = "results/dada2/trackplot.png", width = 7, height = 5)

# summarize percentage of nonchimeric reads as part of merged
nonchim_retained= track %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  select(sample, merged, nonchim) %>%
  group_by(sample) %>%
  mutate(retain = nonchim/merged*100)

x = summary(nonchim_retained$retain) %>% broom::tidy()
write_excel_csv2(x, file = "results/dada2/nochim_sumstat.csv")
# chimera statistics are ok

taxa <- assignTaxonomy(seqtab.nochim, "resources/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "resources/silva_species_assignment_v138.1.fa")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Generate a phylogenetic tree
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences

alignment <- DECIPHER::AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)

write_rds(alignment, file = "intermediate/dada2/decipher_alignment.rds")

# Write DADA2 results to files
saveRDS(taxa, file = "results/dada2/taxa.rds")
saveRDS(fitGTR, file = "results/dada2/fitGTR.rds")

sessionInfo()
