## Analysis of MBVagTrans microbiome data in phyloseq
## Data were provided as fastq files from 16S sequencing
## Primer pair: 341F-805R, Amplicon size 464bp, covering V3-V4 

# Load required packages
library(tidyverse)
library(phyloseq)
library(ggpubr)
library(psadd)
library(DECIPHER)
theme_set(theme_bw())

# Get sample metadata
file_id = readxl::read_xlsx("22010/22010_Analysis/22010_Results.xlsx") %>%
  select(IMGM = "IMGM ID", CID = "Customer-ID", grp = "study group") %>%
  filter(!is.na(IMGM)) %>%
  mutate(grp = factor(grp, labels = c("01_postmenopausal", "02_premenopausal", "03_trans")))
file_id = rbind(file_id, c("22010-PBS01-IKO", 0, "control"))

covar =  read.csv2("22010 Metadaten.csv") %>%
  dplyr::select(CID = "Probenidentifikationsnummer",
                 Age = "Alter.in.Jahren", 
                 GHAT_Length = "Länge.GHAT..in.Monaten", 
                 GNRHA = "GnRh.Analoga.ja.nein", 
                 Virgo, 
                 SexuallyActive = "Sexuell.aktiv", 
                 SexOfPartner = "Biologisches.Geschlecht.Sexualpartner.in",
                 NugentScore = "Nugent.Score", 
                 Testosterone = "Testosteronwert..µg.l.", 
                 Estradiole = "Östradiolwert..ng.l.",
                DurationMenopause = "Zeitraum.seit.Menopause.in.Jahren",
                DurationAmenorrhea = "Zeitraum.seit.letzter.Regelblutung.in.Monaten",
                CycleDaySampling = "Zyklustag.bei.Probenentnahme",
                Duration_GNRH = "Dauer.GnRH.Analoga.Einnahme.in.Monaten") %>%
  dplyr::mutate(GNRHA = factor(GNRHA, labels = c(NA, "yes", "no")),
                Virgo = factor(Virgo, labels = c(NA, "yes", "no")),
                SexuallyActive=ifelse(SexuallyActive=="Ja", "ja", SexuallyActive),
                SexuallyActive = factor(SexuallyActive, labels = c(NA, "yes", "no")),
                SexOfPartner = factor(SexOfPartner, labels = c(NA, "M", "M+F", "F")),
                NugentScore = factor(NugentScore, labels = c("0-3", "4-6", "7-10", "no bacteria")),
                Testosterone = as.numeric(Testosterone),
                Estradiole = as.numeric(Estradiole),
                Age = as.numeric(Age),
                CID = as.character(CID))

# merge sequencing ids and clinical covariates using sample identification number (i.e. customer id CID)
samdat = left_join(file_id, covar, by ="CID")

# save cleaned sample metadata for further downstream use
write_csv2(samdat, file = "intermediate/metadata/mbvagtrans_sampledat_final.csv")

# move sample names into rownames in samdat for importing into phyloseq
samdat2 = column_to_rownames(samdat, "IMGM")

#### Load data from dada2 pipeline ####
seqtab.nochim = readRDS("results/dada2/seqtab_nochim.rds")
taxa = read_rds("results/dada2/taxa.rds")
fitGTR = read_rds("results/dada2/fitGTR.rds")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdat2), 
               tax_table(taxa),
               phy_tree(fitGTR$tree))

ps
set.seed(711)
phy_tree(ps) <- ape::root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
ape::is.rooted(phy_tree(ps))

#### Analysis of PBS Control ####
# export pbs sample into separate phyloseq object and analyse its taxnomic distribution
ps_pbs = prune_samples(sample_names(ps) == "22010-PBS01-IKO", ps)
ps_pbs = prune_taxa(taxa_sums(ps_pbs) > 0, ps_pbs) 
ps_pbs_rel = transform_sample_counts(ps_pbs, function(x) x/sum(x) * 100)
ps_pbs_relmelt = psmelt(ps_pbs_rel) %>%
  mutate(Taxonomy = paste0("p_", Phylum, ";f_", Family, ";g_", Genus)) %>%
  select(Taxonomy, Abundance) %>%
  group_by(Taxonomy) %>%
  summarize(Abundance = sum(Abundance)) %>%
  arrange(-Abundance)
ps_pbs_relmelt
unique(ps_pbs_relmelt$Taxonomy)

write_csv2(ps_pbs_relmelt, file = "results/microbiome/pbs_taxonomic_composition.csv")

# collapse low abundance taxa for plotting
sums <- plyr::ddply(ps_pbs_relmelt, ~Taxonomy, function(x) c(sum=sum(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%
remainder <- sums[sums$sum <= 3,]$Taxonomy

sums[sums$Taxonomy %in% remainder,]$Taxonomy <- 'Other'

sums <- plyr::ddply(sums, ~Taxonomy, function(x) c(sum=sum(x$sum)))
labs <- paste0(round(sums$sum, digits = 1), "%")
pbs_piechart = ggpie(mutate(sums, Abundance = round(sum, digits = 1)),label = labs,
                     lab.font = c(3, "plain", "black"),
                     x="Abundance",
                     fill = "Taxonomy",
                     legend = "right")

ggsave(pbs_piechart, filename = "results/microbiome/pbs_taxonomic_composition_piechart.png")
ggsave(pbs_piechart, filename = "results/microbiome/pbs_taxonomic_composition_piechart.pdf")

#### Clean up ps object ####
ps_samples = prune_samples(sample_names(ps) != "22010-PBS01-IKO", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps_samples))
names(dna) <- taxa_names(ps_samples)
ps_samples <- merge_phyloseq(ps_samples, dna)
taxa_names(ps_samples) <- paste0("ASV", seq(ntaxa(ps_samples)))
ps_samples

# Remove empty taxa from the dataset
ps_samples = prune_taxa(taxa_sums(ps_samples) > 0, ps_samples) # there were no empty taxa

# write the dataset to file, rest of analysis is performed in quarto document
saveRDS(ps_samples, file = "intermediate/ps_samples.rds")

