library(tidyverse)

# load fastq files 
path = "22010/22010_RawData/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# get metadata

# run dada2
library(dada2); packageVersion("dada2")
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path("intermediate", "fq_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("intermediate", "fq_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) 
head(out)

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


plotErrF = plotErrors(errF, nominalQ=TRUE)
plotErrR = plotErrors(errR, nominalQ=TRUE)

saveRDS(errF, file = "intermediate/errF")
saveRDS(errR, file = "intermediate/errR")

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab))) 
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
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

ggplot(track_long, aes(x = step, y = number)) +
  geom_boxplot() +
  geom_jitter()

# summarize percentage of nonchimeric reads as part of merged
nonchim_retained= track %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  select(sample, merged, nonchim) %>%
  group_by(sample) %>%
  mutate(retain = nonchim/merged*100)

summary(nonchim_retained$retain)

# chimera statistics are insuspicious

taxa <- assignTaxonomy(seqtab.nochim, "resources/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "resources/silva_species_assignment_v138.1.fa")

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



#physeq handoff
library(phyloseq); packageVersion("phyloseq")
theme_set(theme_bw())

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps


plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")


