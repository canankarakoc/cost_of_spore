# Spore Costs
# 19 November 2025 - last update
# Author: C. Karakoc
# Reviewed: W. Shoemaker

# Links to the major sources:
# Global expression data from SporeWeb: https://sporeweb.molgenrug.nl/
# Newly synthesized proteins during germination: doi: https://10.1128/mSphere.00463-20 (Swarge et al., 2020)
# Global protein abundance data: https://pax-db.org/
# List of gene categories & annotation: https://subtiwiki.uni-goettingen.de/
# Protein sequence: Uniprot https://www.uniprot.org/taxonomy/224308
# Amino acid/ Nucleotide/ Lipid costs: https://doi.org/10.1073/pnas.1701670114 (Mahmoudabadi et al. 2017)

######################
# Packages & Plotting 
######################
library(stats)
library(tidyverse)       
library(stringr) 
# plotting
library(scales) 
library(ggfortify)
library(ggrepel)
library(ggpubr) 
library(grid)
library(forcats)


mytheme <- theme_bw()+
  theme(axis.ticks.length = unit(.25, "cm"))+
  theme(legend.text = element_text(size=16))+
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18))+
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.border = element_rect(fill=NA, colour = "black", 
size=1))+
theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    linewidth=1)) +
theme(axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
        axis.text.y.right = element_blank(), axis.title.y.right = element_blank())+
   theme(axis.title.x = element_text(margin=margin(10,0,0)),
   axis.title.y = element_text(margin=margin(0,10,0,0)),
   axis.text.x = element_text(margin=margin(10,0,0,0)),
   axis.text.y = element_text(margin=margin(0,10,0,0)))

# Color blind palette
cbpalette <- c("#0072B2", "#D55E00","#009E73", "#CC79A7", "#56B4E9", "#999999", "#F0E442", "#000000")

setwd("~/Documents/GitHub/cost_of_spore")


#######
# Data
#######

# Gene & protein length from SubtiWiki
# Downloaded using wizard of the website 
annotationData <- read.table("./bioaccounting/data/subtiwiki.gene.export.2022-05-11.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Data with synonyms
# Curated from different sources, mainly SubtiWiki synonyms 
nameMap        <- read.table("./bioaccounting/data/nameMap.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Protein abundances from PaxDB
# Downloaded using wizard of the website, data is curated and compiled with manual search using gene synonyms 
protAbun       <- read.table("./bioaccounting/data/protAbunData.csv", sep = ',', dec = ".", header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Protein sequences from Uniprot
# Downloaded using wizard of the website
protSeq        <- read.delim("./bioaccounting/data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.02.08-14.23.07.59.txt")

# Nucleotide & Amino acid costs 
# These datasets are simplified from Mahmoudabadi et al. 2017, Supplementary datasets https://doi.org/10.1073/pnas.1701670114
# See original datasets for comlete breakdown of the costs
aaCosts        <- read.table("./bioaccounting/data/aaCosts_pnas.1701670114.sd01.csv", sep = ",", dec = "." , header = T)
nucCosts       <- read.table("./bioaccounting/data/nucleotideCosts_pnas.1701670114.sd01.csv", sep = ",", dec = "." , header = T)

# Other traits from SubtiWiki lists of lifestyles
otherTraits    <- read.table("./bioaccounting/data/traits.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Germination time course - Swarge et al. 2020 
# This data includes newly synthesized proteins
germination6   <- read.table("./bioaccounting/data/SwargeEtAl_onlyNewProteins.csv", sep = ",", header = T) 
  

                            ################
########################### DATA PREPARATION #####################################
                            ################

#####################
# Protein abundance
#####################

# Protein abundance data does not include locus tags
# It is often easier to merge data sets with locus tags, because genes have a lot of synonyms
# Section below tries to cover as many abundance data as it can, merging with interchangeably used
# genes, synonyms, locus_tags curated by screening

# This is a long pipeline, so sectioned for clarity and readability

# ====================
# Helpers & constants
# ====================
std <- function(x) tolower(gsub("[^a-z0-9]+", "", trimws(x)))
P_total <- 1774445  # proteins per B. subtilis cell (Maass et al. 2011)

# Helper: scientific as "a × 10^b"
sci <- function(x, digits = 1) {
  s <- formatC(x, format = "e", digits = digits)
  sub("e\\+?", " × 10^", s, perl = TRUE)
}

# ==============================================
# SubtiWiki annotation → clean tokens / lengths
# ==============================================
annotation_clean <- annotationData %>%
  mutate(
    locus_tag      = trimws(locus_tag),
    protein_length = suppressWarnings(as.numeric(protein_length)),
    gene_length    = suppressWarnings(as.numeric(gene_length)),
    gene_std       = std(gene)
  )

valid_tags <- annotation_clean %>% distinct(locus_tag)

# Primary gene symbol → locus_tag (priority 2)
primary_map <- annotation_clean %>%
  filter(!is.na(gene_std), gene_std != "") %>%
  distinct(gene_std, locus_tag) %>%
  mutate(score = 2L)

# Ambiguous alias list (priority 1; OK for PaxDB splitting)
alias_cols <- grep("^(geneP|gene[0-9]+|alias[0-9]+)$", names(nameMap), value = TRUE)
alias_map <- if (length(alias_cols) > 0) {
  nameMap %>%
    select(locus_tag, all_of(alias_cols)) %>%
    mutate(across(all_of(alias_cols), ~ na_if(trimws(.), ""))) %>%
    pivot_longer(cols = all_of(alias_cols), names_to = "field", values_to = "alias") %>%
    filter(!is.na(alias)) %>%
    transmute(gene_std = std(alias), locus_tag) %>%
    distinct(gene_std, locus_tag)
} else {
  tibble(gene_std = character(), locus_tag = character())
}

# For *UniProt mapping only*, keep unambiguous versions to avoid 1→many
primary_map_unambig <- primary_map %>%
  add_count(gene_std, name = "n_tags") %>% filter(n_tags == 1L) %>%
  select(gene_std, locus_tag)

alias_map_unambig <- alias_map %>%
  add_count(gene_std, name = "n_tags") %>% filter(n_tags == 1L) %>%
  select(gene_std, locus_tag)

# ======================================================================
# UniProt sequences → locus_tag 
# Priority: BSU token in gene/name (3) > primary_map (2) > alias_map (1)
# ======================================================================
bsu_from_gene <- protSeq %>%
  transmute(protID, locus_tag = str_to_upper(str_extract(gene, "(?i)BSU\\d{5}"))) %>%
  mutate(locus_tag = if_else(nchar(locus_tag) == 8, locus_tag, NA_character_)) %>%
  filter(!is.na(locus_tag)) %>%
  semi_join(valid_tags, by = "locus_tag") %>%
  mutate(score = 3L)

bsu_from_name <- protSeq %>%
  transmute(protID, locus_tag = str_to_upper(str_extract(name, "(?i)BSU\\d{5}"))) %>%
  mutate(locus_tag = if_else(nchar(locus_tag) == 8, locus_tag, NA_character_)) %>%
  filter(!is.na(locus_tag)) %>%
  semi_join(valid_tags, by = "locus_tag") %>%
  mutate(score = 3L)

# Tokenize UniProt "gene" field
protSeq_tok_unique <- protSeq %>%
  mutate(gene_token = strsplit(gene, "[\\s,;/]+")) %>%
  unnest(gene_token) %>%
  mutate(gene_std = std(gene_token)) %>%
  distinct(protID, gene_std)

cand_primary <- protSeq_tok_unique %>%
  inner_join(primary_map_unambig, by = "gene_std") %>%
  mutate(score = 2L) %>%
  select(protID, locus_tag, score) %>%
  distinct(protID, locus_tag, .keep_all = TRUE) %>%
  semi_join(valid_tags, by = "locus_tag")

cand_alias <- protSeq_tok_unique %>%
  inner_join(alias_map_unambig, by = "gene_std") %>%
  mutate(score = 1L) %>%
  select(protID, locus_tag, score) %>%
  distinct(protID, locus_tag, .keep_all = TRUE) %>%
  semi_join(valid_tags, by = "locus_tag")

# Best candidate per UniProt protein
candidates <- bind_rows(bsu_from_gene, bsu_from_name, cand_primary, cand_alias) %>%
  distinct(protID, locus_tag, .keep_all = TRUE)

prot_map <- protSeq %>%
  select(protID, gene, name, sequence) %>%
  inner_join(candidates, by = "protID") %>%
  mutate(seq_len = nchar(sequence)) %>%
  group_by(protID) %>%
  arrange(desc(score), desc(seq_len)) %>%
  slice(1L) %>%
  ungroup() %>%
  select(protID, locus_tag, sequence)

# ====================================================
# PaxDB → locus_tag 
# Split each PaxDB row’s ppm across all matching tags
# ====================================================
# Tokens permitted for mapping: primary symbols, aliases, and BSU IDs
bsu_tokens <- annotation_clean %>%
  transmute(gene_std = std(locus_tag), locus_tag)

tokens_to_tag <- bind_rows(
  primary_map %>% select(gene_std, locus_tag),  # may be many→many
  alias_map   %>% select(gene_std, locus_tag),  # may be many→many
  bsu_tokens                                    # one→one
) %>% distinct()

# Tokenize PaxDB once and de-duplicate tokens per row
protAbun_tok <- protAbun %>%
  mutate(rowid = row_number(),
         gene_token = strsplit(gene, "[\\s,;/]+")) %>%
  unnest(gene_token) %>%
  mutate(gene_std = std(gene_token)) %>%
  distinct(rowid, gene_std, .keep_all = TRUE)

# Split ppm evenly across all matches for each PaxDB row
paxdb_split <- protAbun_tok %>%
  inner_join(tokens_to_tag, by = "gene_std", relationship = "many-to-many") %>%
  add_count(rowid, name = "n_targets") %>%
  mutate(ppm_split = abundance / n_targets) %>%
  group_by(locus_tag) %>%
  summarise(abundance_ppm = sum(ppm_split, na.rm = TRUE), .groups = "drop")

# Coverage diagnostic
sum_paxdb   <- sum(protAbun$abundance, na.rm = TRUE)
sum_matched <- protAbun_tok %>%
  inner_join(tokens_to_tag, by = "gene_std", relationship = "many-to-many") %>%
  distinct(rowid, .keep_all = TRUE) %>%
  summarise(s = sum(abundance, na.rm = TRUE)) %>% pull(s)
coverage_frac <- sum_matched / sum_paxdb
message(sprintf("PaxDB total: %.0f ppm | matched: %.0f ppm (%.1f%%)",
                sum_paxdb, sum_matched, 100*coverage_frac))

# =====================================
# Assemble per-locus_tag protein table 
# =====================================
prot_len_med <- median(annotation_clean$protein_length, na.rm = TRUE)
gene_len_med <- median(annotation_clean$gene_length,    na.rm = TRUE)

protSeqMapped <- prot_map %>%
  left_join(paxdb_split, by = "locus_tag") %>%
  left_join(annotation_clean %>% select(locus_tag, gene, protein_length, gene_length),
            by = "locus_tag") %>%
  mutate(
    protein_length = coalesce(protein_length, nchar(sequence), prot_len_med),
    gene_length    = coalesce(gene_length,    gene_len_med)
  )

# One record per locus_tag (prefer those that actually have mapped ppm, then longest seq)
protSeqTidy <- protSeqMapped %>%
  mutate(has_ppm = !is.na(abundance_ppm), seq_len = nchar(sequence)) %>%
  group_by(locus_tag) %>%
  arrange(desc(has_ppm), desc(seq_len)) %>%
  slice(1L) %>%
  ungroup()

# ============================
# Amino-acid cost per protein
# ============================
alphabet <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

seqCount <- protSeqTidy %>%
  select(locus_tag, sequence) %>%
  rowwise() %>%
  reframe(
    locus_tag = locus_tag,
    symbol    = alphabet,
    aac       = str_count(sequence, alphabet),
    .groups   = "drop_last"
  ) %>%
  ungroup() %>%
  left_join(aaCosts, by = "symbol") %>%
  summarise(
    aa_opportunitySum = sum(aac * opportunity_costs, na.rm = TRUE),
    aa_directSum      = sum(aac * direct_costs,      na.rm = TRUE),
    .by = locus_tag
  )

# ===================================================================================
#  Final tidy tables for downstream
#    - _base: raw mapped ppm, no imputation
#    - main: adds tiny floor + converts to absolute copies per cell, then back to ppm
# ===================================================================================
protSeqTidyAbun_base <- protSeqTidy %>%
  left_join(seqCount, by = "locus_tag") %>%
  mutate(abundance_ppm_raw = suppressWarnings(as.numeric(abundance_ppm))) %>%
  select(protID, locus_tag, gene, abundance_ppm_raw, protein_length, gene_length,
         aa_opportunitySum, aa_directSum)

# Use a very small ppm floor only where missing, then scale to P_total proteins/cell
floor_ppm <- 0.1  # tweakable; negligible vs top abundances
protSeqTidyAbun <- protSeqTidyAbun_base %>%
  mutate(
    abundance_ppm_used  = coalesce(abundance_ppm_raw, floor_ppm),
    # keep absolute counts sensible: distribute P_total by relative ppm weights
    proteins_per_cell   = (abundance_ppm_used / sum(abundance_ppm_used, na.rm = TRUE)) * P_total,
    abundance_ppm_final = proteins_per_cell / P_total * 1e6,   # closed composition for ppm outputs
    # for compatibility with older code that expects a column named `abundance` in ppm:
    abundance           = abundance_ppm_final
  )

#################################################
# Expression data merged with protein abundances 
#################################################

# Expression data from SporeWeb

files  <- list.files(path = "./bioaccounting/data/SporeWebHeatmaps/" , pattern = "*.csv", full.names = T)
exp_files <-  list()
for (i in 1:length(files)){
  exp_files[[i]] <- read.table(files[i], header = T, sep = ",", dec = ".")
}
merged_exp                   <- bind_rows(exp_files, .id = 'sourceID')
merged_exp$sourceID          <- as.factor(merged_exp$sourceID)
levels(merged_exp$sourceID ) <- c("1.vegetative", "2.starvation", "3.onset", "4.commitment", "5.engulfment")

# Long formate of expression 
expressionLong <- merged_exp %>%
  select(gene, locus_tag, regulators, sourceID, starts_with("t")) %>%
  pivot_longer(cols = starts_with("t"),
               names_to = "time", values_to = "expression") %>%
  group_by(time, sourceID, regulators, gene, locus_tag) %>%
  summarise(meanexp = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  filter(meanexp > 0) %>%
  mutate(time_h = sub("^t", "", time))

# Spore gene set (unique locus_tags that appear in the time series) 
mergedExpData_time <- expressionLong %>%
  left_join(protSeqTidyAbun, by = "locus_tag")   # keeps all time points

gene_len_med <- median(annotation_clean$gene_length, na.rm = TRUE)

spore_once <- mergedExpData_time %>%
  filter(!is.na(locus_tag)) %>%
  distinct(locus_tag, .keep_all = TRUE) %>%      # one row per gene
  transmute(
    locus_tag,
    gene_length = suppressWarnings(as.numeric(gene_length)) 
  ) %>%
  mutate(gene_length = coalesce(gene_length, gene_len_med))

#################################################
# Germination data merged with protein abundances 
#################################################

# Interval table from Swarge et al.2020 (differences between successive time points)
# Note in the paper after 210 min is considered as vegetative growth
# I am excluding the last point not to double the costs. 

# Keep only up to 3.5 h
keep_intervals <- c("H0.25","H0.5","H1","H1.5","H2.5","H3.5")

germ_intervals <- germination6 %>%
  rename(protID = ProtID) %>%
  transmute(
    protID,
    H0.25 = T15  - T0,
    H0.5  = T30  - T15,
    H1    = T60  - T30,
    H1.5  = T90  - T60,
    H2.5  = T150 - T90,
    H3.5  = T210 - T150
    # H5.5 (T330-T210) intentionally omitted
  ) %>%
  pivot_longer(starts_with("H"),
               names_to = "interval", values_to = "score") %>%
  mutate(score = pmax(score, 0)) %>%   # clamp any negatives to 0
  filter(interval %in% keep_intervals, score > 0)

# Interval metadata (durations and midpoints for plotting)
interval_meta <- tibble(
  interval = keep_intervals,
  dur_h    = c(0.25, 0.25, 0.50, 0.50, 1.0, 1.0),
  mid_h    = c(0.125, 0.375, 0.75, 1.25, 2.0, 3.0)
)

germ_w <- germ_intervals %>%
  left_join(interval_meta, by = "interval") %>%
  left_join(
    protSeqTidyAbun %>%
      select(protID, locus_tag, proteins_per_cell,
             gene_length, protein_length, aa_directSum, aa_opportunitySum),
    by = "protID"
  ) %>%
  group_by(protID) %>%
  mutate(
    w = score / sum(score, na.rm = TRUE),  # renormalized over kept bins
    w = if_else(is.finite(w) & w > 0, w, 0)
  ) %>%
  ungroup()

# One row per germination/outgrowth gene, with safe gene_length
gene_len_med <- if (exists("gene_len_med")) gene_len_med else median(annotation_clean$gene_length, na.rm = TRUE)

germ_once <- germ_w %>%
  dplyr::filter(w > 0, !is.na(locus_tag)) %>%
  dplyr::distinct(locus_tag, .keep_all = TRUE) %>%
  dplyr::transmute(
    locus_tag,
    gene_length = suppressWarnings(as.numeric(gene_length))
  ) %>%
  dplyr::mutate(gene_length = dplyr::coalesce(gene_length, gene_len_med))


########################################
# Whole genome and whole membrane costs
# Necessary for downstream analysis
# Used in multiple sections
########################################

#______________________________________________________________
# Replication costs (Whole genome) 
# Paid during spore formation
# Manuscript: 
# Results:"Energetics of spore formation", Figure 1 (pie chart)
# Results: "Head-to-head comparison", Figure 3 (genome bar)
# Methods: "Bioenergetic accounting: definitions and assumptions", 
# "Replication costs" 
# Equations: 5-6
#______________________________________________________________

# DNA unwinding 1 ATP per base pair # Lee and Yang, 2006 (https://doi.org/10.1016/j.cell.2006.10.049)
# Primer of Okazaki fragments 0.32 ATP #Lynch and Marinov, 2015
# Ligation is negligible   

# Constants 
genome.size <- 4215606 #NCBI Reference Sequence: NC_000964.3
# (Lg notation in the equations) 

# Nucleotide costs 
dna_nucleotide <- nucCosts %>% filter(molecule == "DNA")
rna_nucleotide <- nucCosts %>% filter(molecule == "RNA")

dna_nucleotide_cost_opp <- mean(dna_nucleotide$opportunity) # (c_o)
dna_nucleotide_cost_dir <- mean(dna_nucleotide$direct)  # (c_s + c_p), includes polymerization costs 2ATPs

rna_nucleotide_cost_opp <- mean(rna_nucleotide$opportunity)
rna_nucleotide_cost_dir <- mean(rna_nucleotide$direct) # includes polymerization costs 2ATPs

# Other constant components
helicase_cost <- 1 # (c_hel)
primers_cost  <- 0.32 # lagging-strand primers (c_prim)

# Opportunity
genome_opp <- 2 * genome.size * dna_nucleotide_cost_opp  # 286661208 PO

# Direct 
genome_dir <- 2 * genome.size * (dna_nucleotide_cost_dir) + (helicase_cost * genome.size) + (primers_cost * genome.size)  # 114116454 PD

# Total
genome_tot <- genome_opp + genome_dir # 400777662 PT

#_________________________________________________________
# Total Membrane lipid costs
# Manuscript: 
# Results: "Head to head comparison", Figure 3 (membrane bar)
# Methods: "Bioenergetic accounting: definitions and assumptions, 
# Membrane synthesis and remodelling"
# Equations: 12-14
#_________________________________________________________

# Number of lipid molecules = Cellular membrane areas/head-group areas of membrane lipid molecules
# Head group area is a1 = 0.65 nm2 (Nagle and Tristram-Nagle 2000; Petrache et al. 2000; Kucerka et al. 2011).

# Thickness of the bilayer (h):
# The thickness of a bilayer is approximately twice the radius of the head-group area, which 0.5 nm in all cases,
# plus the total length of the internal hydrophobic tail domains (Lewis and Engelman 1983; Mitra et al. 2004), 
# generally 3.0 nm, so total is 4 nm.

# Bacillus average length (a) and width (b) Barak et al. 2018. 

# constants
a1 <- 0.65
h  <- 4
L  <- 2.5e3
D  <- 1.0e3
a  <- L/2
b  <- D/2

stopifnot(h > 0, b > h, L > 2*h)

protein_area_frac <- 0.5
lipid_opp <- 212 
lipid_dir <- 18

# spherocylinder surfaces
A_outer <- 4*pi*a*b               # = 2*pi*b*L
A_inner <- 4*pi*(a - h)*(b - h)   # = 2*pi*(b-h)*(L-2*h)

N_lipid_mem <- (A_outer + A_inner) / a1 * (1 - protein_area_frac)

membraneOpp <- N_lipid_mem * lipid_opp   # 2547294111 PO
membraneDir <- N_lipid_mem * lipid_dir   # 216279689 PD
membraneTot <- membraneOpp + membraneDir # 2763573799 PT


                            ################
########################### SPORE FORMATION  #####################################
                            ################

#______________________________________________________________________
# Replication costs of expressed genes during spore formation
# Relative costs compared to to the whole genome
# Manuscript: 
# Results: "Energetics of spore formation", first paragraph
# Same equations as the whole genome synthesis Eq.s 5-6
#______________________________________________________________________

spore_related_genes_length <- sum(spore_once$gene_length, na.rm = TRUE)
spore_related_genes_number <- nrow(spore_once)
total_bacillus_genes       <- 4237 # NCBI Genes annotated on Bacillus subtilis subsp. subtilis str. 168 ASM904v1 (GCF_000009045.1)

rel_genome <- spore_related_genes_length / genome.size   # length fraction %18
rel_genes  <- spore_related_genes_number / total_bacillus_genes  # gene-count fraction %21

# Costs of spore formation related genes only  
sporeRep <- spore_once %>%
  mutate(
    opportunity = 2 * gene_length * dna_nucleotide_cost_opp,
    direct      = 2 * gene_length * dna_nucleotide_cost_dir +
      gene_length * helicase_cost +
      gene_length * primers_cost
  )

sporeRepSum <- summarise(sporeRep,
                         sumOpp = sum(opportunity, na.rm = TRUE),
                         sumDir = sum(direct,      na.rm = TRUE))

sporeRepTotal <- sporeRepSum$sumOpp + sporeRepSum$sumDir # 51307700 PO, 20424992 PD, 71732692 PT
sporeRepTotal_frac_of_genome <- sporeRepTotal / genome_tot # 18%

#_______________________________________________________________
# Spore Transcription
# Manuscript: 
# Results: "Energetics of spore formation", first paragraph, 
# Figure 1 (PD, PO bars)
# Methods: "Bioenergetic accounting: definitions and assumptions, 
# "Transcription costs" 
# Equations: 7-9
#_______________________________________________________________

# 1 mRNA can yield to 100-1000 proteins (Cell biology by the numbers)
# Protein abundance is reported as parts per million (ppm). 
# Number of proteins was reported experimentally as average 1774445 in Bacillus (Maass et.al. 2011, DOI: 10.1021/ac1031836)
# Average sporulation time is 8hours, median mRNA degradation rate of Bacillus is 12 (mRNAs degraded every 5 minutes) per hour 
# Hambraeus et a., 2003 (DOI: 10.1007/s00438-003-0883-6), assuming nucleotides are well recycled and it only affects polymerization costs
# This rate is not used numerically here, because we do not have mRNA copy numbers, we approximate from protein demand (ppm), 
# and expression weights, which yields identical totals without double-counting
# I account synthesis costs separately for direct costs, so I can consider repolymerization costs to keep up with mRNA copies

poly_costs <- 2  # per nt (c_pr)
rna_nucleotide_cost_synthesis <- rna_nucleotide_cost_dir - poly_costs  # (c_sr)

yield   <- 100 # (Y in the manuscript)
ppm_to_count <- function(ppm) (ppm / 1e6) * P_total

# tiny extra direct cost per nt for gyrase/supercoil removal (tunable; 0.1 is an upper bound)
c_supercoil <- 0.1 # PD per nt (c_sc)

# weights by hour (use expression as weights; safe fallback if zero/NA)
expr_w <- mergedExpData_time %>%
  select(locus_tag, time_h, meanexp, proteins_per_cell, gene_length) %>%
  group_by(locus_tag) %>%
  mutate(
    w = meanexp / sum(meanexp, na.rm = TRUE),
    w = if_else(is.finite(w) & !is.na(w), w, 1 / n())
  ) %>%
  ungroup()

# first appearance of each gene 
first_appearance <- mergedExpData_time %>%
  group_by(locus_tag) %>%
  slice_min(as.numeric(time_h), with_ties = FALSE) %>%
  ungroup() %>%
  transmute(locus_tag, time_h_first = time_h)

# As described in the manuscript: 
# (i) One-time synthesis + opportunity at first appearance (C_TD (synth) and C_TO)
one_time <- expr_w %>%
  distinct(locus_tag, proteins_per_cell, gene_length) %>%
  mutate(
    transcripts_total = proteins_per_cell / yield,
    nt_total          = gene_length * transcripts_total
  ) %>%
  left_join(first_appearance, by = "locus_tag") %>%
  group_by(time_h_first) %>%
  summarise(
    dir_synth_once = sum(nt_total * rna_nucleotide_cost_synthesis, na.rm = TRUE), # (c_sr)
    opp_once       = sum(nt_total * rna_nucleotide_cost_opp,       na.rm = TRUE), # (c_or)
    .groups = "drop"
  )

# (ii) Hourly polymerization from first appearance onward (C_TD (t))
expr_w_post <- expr_w %>%
  left_join(first_appearance, by = "locus_tag") %>%
  mutate(
    time_num  = as.numeric(time_h),
    first_num = as.numeric(time_h_first),
    w = if_else(time_num < first_num, 0, w)
  ) %>%
  group_by(locus_tag) %>%
  mutate(w = if (sum(w, na.rm = TRUE) > 0) w / sum(w, na.rm = TRUE) else 0) %>%
  ungroup()

per_hour_terms <- expr_w_post %>%
  mutate(
    transcripts_total = proteins_per_cell / yield,
    transcripts_hour  = transcripts_total * w,
    nt_hour           = gene_length * transcripts_hour,
    dir_poly_hour     = nt_hour * poly_costs,    # 2 PD/nt # (c_pr)
    dir_super_hour    = nt_hour * c_supercoil    # small extra PD/nt # (c_sc)
  ) %>%
  group_by(time_h) %>%
  summarise(
    dir_poly_hour  = sum(dir_poly_hour,  na.rm = TRUE),
    dir_super_hour = sum(dir_super_hour, na.rm = TRUE),
    .groups = "drop"
  )

# Combine to per-hour totals
transcription_hourly <- per_hour_terms %>%
  full_join(one_time, by = c("time_h" = "time_h_first")) %>%
  replace_na(list(dir_synth_once = 0, opp_once = 0,
                         dir_poly_hour = 0,  dir_super_hour = 0)) %>%
  mutate(
    direct      = dir_synth_once + dir_poly_hour + dir_super_hour,
    opportunity = opp_once,
    total       = direct + opportunity
  ) %>%
  arrange(as.numeric(time_h))

total_direct <- sum(transcription_hourly$direct)      # 21756135 PD
total_opp    <- sum(transcription_hourly$opportunity) # 56637872 PO
total_all    <- sum(transcription_hourly$total)       # 78394007 PT

#_______________________________________________________________
# Spore Translation
# Manuscript: 
# Results: "Energetics of spore formation", first paragraph, 
# Figure 1 (PD, PO bars)
# Methods: "Bioenergetic accounting: definitions and assumptions, 
# "Translational costs" 
# Equations: 10-11
#_______________________________________________________________

# Base table: only genes that appear in the time series (spore set)
genes_once_aaAbun <- protSeqTidyAbun %>%
  select(locus_tag, gene,
                proteins_per_cell,     # absolute copies per cell (already scaled)
                protein_length,        # numeric (filled earlier)
                aa_directSum,          # PD per protein (from sequence) (LC_D)
                aa_opportunitySum) %>% # PO per protein (from sequence) (LC_O)
  inner_join(first_appearance, by = "locus_tag")   # drop non-spore genes

# ~11 of the “spore genes” don’t have a corresponding protein row in protSeqTidyAbun. Common reasons:
# the locus_tag is non–protein-coding (ncRNA, asRNA, riboswitch, etc.),
# it’s a protein but didn’t map to UniProt (no reviewed entry / naming mismatch),
# or a stray locus_tag variant.

# Imputation strategy for missing AA sums:
# i) prefer length * median(per-aa cost) across proteins
# ii) if length is NA (should be rare), fall back to protein-level medians
dir_per_aa_med <- median(genes_once_aaAbun$aa_directSum      /
                           genes_once_aaAbun$protein_length, na.rm = TRUE)
opp_per_aa_med <- median(genes_once_aaAbun$aa_opportunitySum /
                           genes_once_aaAbun$protein_length, na.rm = TRUE)

med_aa_dir <- median(genes_once_aaAbun$aa_directSum,      na.rm = TRUE)
med_aa_opp <- median(genes_once_aaAbun$aa_opportunitySum, na.rm = TRUE)

translation_by_gene <- genes_once_aaAbun %>%
  mutate(
    aa_direct_imputed = coalesce(
      aa_directSum,
      protein_length * dir_per_aa_med,
      med_aa_dir
    ),
    aa_opp_imputed = coalesce(
      aa_opportunitySum,
      protein_length * opp_per_aa_med,
      med_aa_opp
    ),
    direct      = proteins_per_cell * aa_direct_imputed, # (C_Pr_D)
    opportunity = proteins_per_cell * aa_opp_imputed,    # (C_Pr_O)
    total       = direct + opportunity
  )

## Translation extras (direct P_D only)
## init + term ≈ 1 GTP each → ~2 P_D per protein

# toggle / params
add_init_term <- TRUE
init_GTP <- 1   # IF2
term_GTP <- 1   # RF3
init_term_PD_per_protein <- (init_GTP + term_GTP) * 1  # 1 GTP ~ 1 P_D

# Per-gene extras applied
translation_by_gene_extras <- translation_by_gene %>%
  mutate(
    # base direct from earlier step
    direct_base = direct,
    
    # extras (all P_D)
    pd_init_term = if (isTRUE(add_init_term)) proteins_per_cell * init_term_PD_per_protein else 0,
    
    # updated direct/total
    direct_adj = direct_base + pd_init_term,
    total_adj  = direct_adj + opportunity
  )

# Time series with extras applied (booked at first appearance)
translation_by_time_extras <- translation_by_gene_extras %>%
 group_by(time_h_first) %>%
  summarise(
    direct      = sum(direct_adj,  na.rm = TRUE),
    opportunity = sum(opportunity, na.rm = TRUE),
    total       = sum(total_adj,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(time_h = time_h_first) %>%
  arrange(as.numeric(time_h))

# Totals
translation_totals_with_extras <- translation_by_time_extras %>%
  summarise(
    total_direct      = sum(direct),      # 320262066 PD
    total_opportunity = sum(opportunity), # 1308452259 PO
    total_all         = sum(total)        # 1628714325 PT
  )

#___________________________________________________________
# Septum synthesis
# Manuscript:
# Results: "Energetics of spore formation", first paragraph, 
# Figure 1 (pie chart)
# Methods: "Membrane synthesis and remodelling"
# Equations: 15
#___________________________________________________________

# septum as a flat circular bilayer 
# apply same 50% protein-occlusion discount to the septal bilayer
# Septum as a flat circular bilayer (diameter ~ D = 1 µm, radius b = D/2)
A_sept_bilayer <- 2 * pi * b^2  # both leaflets, same planar area
N_lipid_sept   <- (A_sept_bilayer / a1) * (1 - protein_area_frac)

septumOpp <- N_lipid_sept * lipid_opp # 256160632 PO
septumDir <- N_lipid_sept * lipid_dir # 21749488 PD
septumTot <- septumOpp + septumDir    # 277910119 PT

#______________________________
# Total costs, plots, pie, bars
# Following will yield Figure 1
#______________________________

## Hourly totals (transcription + translation) 
hour_grid <- tibble(time_h = as.character(1:8))

hourly <- hour_grid %>%
  left_join(
    transcription_hourly %>% select(time_h,
                                           trscr_direct = direct,
                                           trscr_opp    = opportunity),
    by = "time_h"
  ) %>%
  left_join(
    translation_by_time_extras %>% select(time_h,
                                          transl_direct = direct,
                                          transl_opp    = opportunity),
    by = "time_h"
  ) %>%
  replace_na(list(trscr_direct = 0, trscr_opp = 0,
                         transl_direct = 0, transl_opp = 0)) %>%
  mutate(
    direct      = trscr_direct + transl_direct,
    opportunity = trscr_opp    + transl_opp,
    total       = direct + opportunity
  ) %>%
  arrange(as.numeric(time_h))

## Pie components 
cost_rep_all    <- genome_opp + genome_dir     # genome_tot
cost_transcript <- sum(total_direct + total_opp, na.rm = TRUE)
cost_translation<- sum(translation_totals_with_extras$total_direct + translation_totals_with_extras$total_opportunity, na.rm = TRUE)

# Membrane 
cost_membrane <- septumTot # (septumOpp + septumDir)

# Ultimate total with estimated membrane costs 
all_pie_costs <- cost_rep_all + cost_transcript + cost_translation + cost_membrane #2385616310 PT

pieData <- tibble(
  pieCost     = c("replication", "transcription", "translation", "membrane"),
  proportion  = 100 * c(cost_rep_all, cost_transcript, cost_translation, cost_membrane) / all_pie_costs
) %>%
  mutate(labels = paste0(pieCost, " ", round(proportion, 1), "%"))

# Plot (uses your mytheme + ggrepel)
pieSpore <- ggplot(pieData, aes(x = "", y = proportion, fill = pieCost))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  mytheme +
  scale_fill_manual(values = c( "#CC79A7", "#009E73", "#D55E00", "#0072B2")) +
  geom_text_repel(aes(label = labels), size = 4.5, show.legend = FALSE) +
  theme_void() +
  theme(legend.position = "none")

#ggsave("bioaccounting/figures/figure1_pieSpore.pdf", pieSpore, height = 5, width = 6)

## Opportunity/Direct ratio 
# Manuscript "Energetics of spore formation", first paragraph
opportunity_all <- sum(genome_opp + total_opp +  translation_totals_with_extras$total_opportunity + septumOpp)
direct_all      <- sum(genome_dir + total_direct +  translation_totals_with_extras$total_direct + septumDir)
 
opportunity_all/all_pie_costs #80% 
direct_all/all_pie_costs #20% 

## Costs of first hour, relative to total costs, including replication
# Manuscript "Energetics of spore formation", second paragraph
first_hour <- (as.numeric(hourly %>% filter(time_h == 1) %>% select(total))) + cost_rep_all
first_hour/all_pie_costs #67%

## Time-series (stacked bars + asymptotic fit) 
sporulationCosts <- rbind(
  transmute(hourly, time = as.numeric(time_h), type = "opportunity", costs = opportunity),
  transmute(hourly, time = as.numeric(time_h), type = "direct",      costs = direct)
)

sporulationCosts_lay1 <- hourly %>%
  transmute(time = as.numeric(time_h), sum = total)

# Nonlinear asymptotic fit to total
asym_data <- transmute(hourly, time = as.numeric(time_h), total = total)

# Fit
fit <- nls(total ~ SSasymp(time, yf, y0, log_alpha), data = asym_data)

# Parameters
pars   <- coef(fit)
lambda <- exp(pars[["log_alpha"]])   # SSasymp uses log_alpha

# Predictions for curve
tt   <- seq(1, 8, by = 0.1)
pred <- predict(fit, newdata = data.frame(time = tt))
preddata <- data.frame(tt = tt, pred = pred)

## R^2 on original scale
y    <- asym_data$total
yhat <- predict(fit)
RSS  <- sum((y - yhat)^2)
TSS  <- sum((y - mean(y))^2)
R2   <- 1 - RSS/TSS

# Annotations
lab_lambda <- paste0("\u03BB = ", round(lambda, 3))   # λ
lab_R2     <- paste0("R^2 = ", round(R2, 3))

my_lab <- c(expression(P['D']), expression(P['O']), expression(P['T']))

f1 <- ggplot() +
  geom_vline(xintercept = 2, linetype = "dashed") +
  geom_bar(data = sporulationCosts_lay1,
                    aes(x = time, y = sum),
                    stat = "identity", color = "grey90", fill = "grey75", alpha = 0.5) +
  geom_bar(data = sporulationCosts,
                    aes(x = time, fill = type, y = costs),
                    stat = "identity", color = "grey25", position = position_dodge(width = 1)) +
  ylab("ATP molecules") +
  xlab("Time (h)") +
  geom_line(data = preddata, aes(x = tt, y = pred), linewidth = 1) +
  mytheme +
  scale_y_continuous(breaks = c(2e8, 4e8, 6e8, 8e8, 1e9),
                              labels = c(2, 4, 6, 8, 10),
                              sec.axis = dup_axis()) +
  scale_x_continuous(sec.axis = dup_axis()) +
  annotate("text", x = -0.7, y = 1.2e9, label = "(x10^8)", parse = TRUE, size = 18/.pt) +
  coord_cartesian(xlim = c(0.5, 8.5), clip = "off") +
  annotate("text", x = 8, y = 2e8,
                    label = paste0("lambda==", round(lambda, 3)),
                    hjust = "right", size = 6, fontface = 'italic', parse = TRUE) +
  annotate("text", x = 8, y = 1.2e8,
                    label = paste0("R^2==", round(R2, 3)),
                    hjust = "right", size = 6, fontface = 'italic', parse = TRUE) +
  theme(legend.position = c(0.32, 0.83), legend.title = element_blank()) +
  scale_fill_manual(values = c("#D55E00", "#0072B2"),
                    labels = c(expression(P[D]), expression(P[O])))

# ggsave("bioaccounting/figures/figure1_sporeCostsTime.pdf", f1, height = 5, width = 6)


## ------------------------------------------
## Practical summary for the manuscript text
# Checks if the code and manuscript match
## ------------------------------------------

# Component PTs (for process shares)
cost_rep_all    <- genome_tot
cost_transcript <- total_direct + total_opp
cost_translation<- translation_totals_with_extras$total_direct +
  translation_totals_with_extras$total_opportunity
cost_membrane   <- septumTot

spore_total_PT  <- all_pie_costs
spore_PD_total  <- direct_all
spore_PO_total  <- opportunity_all
spore_PD_pct    <- 100 * spore_PD_total / spore_total_PT
spore_PO_pct    <- 100 * spore_PO_total / spore_total_PT

process_shares <- tibble::tibble(
  process = c("translation","genome replication","membrane (septum)","transcription"),
  PT      = c(cost_translation, cost_rep_all, cost_membrane, cost_transcript)
) %>%
  dplyr::mutate(share_pct = 100 * PT / spore_total_PT)

# sanity: shares sum to ~100%
stopifnot(abs(sum(process_shares$share_pct) - 100) < 1e-6)

# Table (PT, PD, PO) for the four bins + grand total
ms_spore_summary <- tibble::tribble(
  ~item,                   ~PT,              ~PD,                                                     ~PO,                               ~note,
  "Translation",           cost_translation,  translation_totals_with_extras$total_direct,            translation_totals_with_extras$total_opportunity, "all spore proteins",
  "Transcription",         cost_transcript,   total_direct,                                           total_opp,                         "spore genes only",
  "Genome replication",    cost_rep_all,      genome_dir,                                             genome_opp,                        "whole chromosome duplicated",
  "Septum lipids",         cost_membrane,     septumDir,                                              septumOpp,                         "1 µm disk, 50% occlusion",
  "Total (sum above)",     spore_total_PT,    spore_PD_total,                                         spore_PO_total,                    ""
) %>%
  dplyr::mutate(
    PT_str = sci(PT, 1),
    PD_str = sci(PD, 1),
    PO_str = sci(PO, 1)
  ) %>%
  dplyr::select(item, PT_str, PD_str, PO_str, note)

# Headline helpers for prose
headline <- list(
  total_PT_str      = sci(spore_total_PT, 1),
  PD_pct            = round(spore_PD_pct, 1),
  PO_pct            = round(spore_PO_pct, 1),
  share_translation = round(process_shares$share_pct[process_shares$process=="translation"], 1),
  share_replication = round(process_shares$share_pct[process_shares$process=="genome replication"], 1),
  share_membrane    = round(process_shares$share_pct[process_shares$process=="membrane (septum)"], 1),
  share_transcript  = round(process_shares$share_pct[process_shares$process=="transcription"], 1),
  genome_PT_str     = sci(genome_tot, 1),
  spore_gene_frac   = round(100*rel_genes, 1),
  spore_rep_frac    = round(100*sporeRepTotal_frac_of_genome, 1),
  septum_pct_of_PT  = round(100*cost_membrane/spore_total_PT, 2)
)

# Print block
cat("\n=== Spore formation summary (PT, PD, PO) ===\n")
print(ms_spore_summary, n = Inf)

cat("\n=== Shares of total spore PT (%) ===\n")
print(process_shares %>% dplyr::mutate(share_pct = round(share_pct, 1)) %>% dplyr::select(process, share_pct))

cat(sprintf(
  "\nText helpers:\nTotal spore cost ≈ %s ATP; partition ~%s%% P_O / %s%% P_D.\nBreakdown: TL %s%%, Rep %s%%, Membrane %s%%, TX %s%%.\nGenome duplication ≈ %s ATP.\nSporulation genes ≈ %s%% of genome and account for ≈ %s%% of replication cost.\nSeptum ≈ %s%% of P_T.\n",
  headline$total_PT_str, headline$PO_pct, headline$PD_pct,
  headline$share_translation, headline$share_replication,
  headline$share_membrane, headline$share_transcript,
  headline$genome_PT_str, headline$spore_gene_frac,
  headline$spore_rep_frac, headline$septum_pct_of_PT
))

###########################################
# Other numbers reported in the manuscript:
###########################################

#________________________________________________________________________________
# What's the percentage of mother cell and forespore contribution to the costs?
# Manuscript:
# Results: "Energetics of spore formation", second paragraph
#________________________________________________________________________________

# Map regulators -> compartment 
reg_class <- mergedExpData_time %>%
  select(locus_tag, regulators) %>%
  mutate(reg = tolower(regulators)) %>%
  # presence of sigma factors (sign +/− doesn't matter for tagging)
  transmute(
    locus_tag,
    hasE = stringr::str_detect(reg, "sige"),
    hasK = stringr::str_detect(reg, "sigk"),
    hasF = stringr::str_detect(reg, "sigf"),
    hasG = stringr::str_detect(reg, "sigg")
  ) %>%
  group_by(locus_tag) %>%
  summarise(
    hasE = any(hasE, na.rm = TRUE),
    hasK = any(hasK, na.rm = TRUE),
    hasF = any(hasF, na.rm = TRUE),
    hasG = any(hasG, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mother = hasE | hasK,          # σE/σK => mother cell
    spore  = hasF | hasG,          # σF/σG => forespore
    compartment = case_when(
      mother & !spore ~ "mother",
      spore  & !mother ~ "forespore",
      mother &  spore  ~ "both",
      TRUE             ~ "unknown"
    )
  ) 

# Per-gene totals for TX and TL (PD/PO/PT)
# Translation per gene
tl_by_gene <- translation_by_gene_extras %>%
  select(locus_tag,
                TL_PD = direct_adj,
                TL_PO = opportunity) %>%
  mutate(TL_PT = TL_PD + TL_PO)

# Transcription per gene: one-time synthesis + total polymerization (+tiny supercoil)
tx_by_gene <- expr_w %>%                         # from your earlier pipeline
  distinct(locus_tag, proteins_per_cell, gene_length) %>%
  mutate(
    transcripts_total = proteins_per_cell / yield,
    nt_total = gene_length * transcripts_total,
    TX_PO = nt_total * rna_nucleotide_cost_opp,
    TX_PD = nt_total * (rna_nucleotide_cost_synthesis + poly_costs + c_supercoil),
    TX_PT = TX_PD + TX_PO
  ) %>%
  select(locus_tag, TX_PD, TX_PO, TX_PT)

prog_gene_costs <- tx_by_gene %>%
  left_join(tl_by_gene, by = "locus_tag") %>%
  replace_na(list(TL_PD = 0, TL_PO = 0, TL_PT = 0)) %>%
  mutate(
    PD = TX_PD + TL_PD,
    PO = TX_PO + TL_PO,
    PT = PD + PO
  ) %>%
  select(locus_tag, PD, PO, PT)


# Mother vs forespore totals
# exclude "both" and "unknown"
by_comp_simple <- prog_gene_costs %>%
  inner_join(reg_class, by = "locus_tag") %>%
  filter(compartment %in% c("mother", "forespore")) %>%
  group_by(compartment) %>%
  summarise(PT = sum(PT, na.rm = TRUE),
                   PD = sum(PD, na.rm = TRUE),
                   PO = sum(PO, na.rm = TRUE),
                   .groups = "drop")

# split "both" genes 50/50 between mother & forespore
both_half <- prog_gene_costs %>%
  inner_join(reg_class, by = "locus_tag") %>%
  filter(compartment == "both") %>%
  mutate(weight = 0.5) %>%
  select(locus_tag, PD, PO, PT, weight) %>%
  {bind_rows(mutate(., compartment = "mother"),
                     mutate(., compartment = "forespore")) }

single_full <- prog_gene_costs %>%
  inner_join(reg_class, by = "locus_tag") %>%
  filter(compartment %in% c("mother", "forespore")) %>%
  mutate(weight = 1)

by_comp_split <- bind_rows(single_full, both_half) %>%
  group_by(compartment) %>%
  summarise(PT = sum(PT * weight, na.rm = TRUE),
                   PD = sum(PD * weight, na.rm = TRUE),
                   PO = sum(PO * weight, na.rm = TRUE),
                   .groups = "drop")

# Report for mannuscript
report_shares <- function(df) {
  tot <- sum(df$PT)
  mom <- df %>% filter(compartment == "mother") %>% dplyr::pull(PT)
  fore<- df %>% filter(compartment == "forespore") %>% dplyr::pull(PT)
  cat(sprintf("Mother: %.1f%%   Forespore: %.1f%% (of TX+TL PT)\n",
              100*mom/tot, 100*fore/tot))
  df
}

cat("Excluding ambiguous genes:\n")
print(report_shares(by_comp_simple))

cat("\nSplitting 'both' genes 50/50:\n")
print(report_shares(by_comp_split))

# TX+TL split
tx_tl_comp <- by_comp_split  # columns: compartment, PD, PO, PT

# Fixed costs you to assign to the mother cell
#    (whole-genome replication + septum lipids during sporulation)
fixed_to_mother <- tibble::tibble(
  compartment = "mother",
  PD = genome_dir + septumDir,
  PO = genome_opp + septumOpp,
  PT = (genome_dir + genome_opp) + (septumDir + septumOpp)
)

# Combine and summarise
comp_with_fixed <- dplyr::bind_rows(
  dplyr::mutate(tx_tl_comp, source = "TX+TL"),
  dplyr::mutate(fixed_to_mother, source = "replication+septum")
) %>%
  dplyr::group_by(compartment) %>%
  dplyr::summarise(PD = sum(PD, na.rm = TRUE),
                   PO = sum(PO, na.rm = TRUE),
                   PT = sum(PT, na.rm = TRUE),
                   .groups = "drop")

# Report
share <- function(df) {
  tot <- sum(df$PT, na.rm = TRUE)
  df %>%
    dplyr::transmute(compartment,
                     PT       = PT,
                     PD       = PD,
                     PO       = PO,
                     share_PT = 100 * PT / tot)
}

cat("\n— TX+TL only —\n")
print(share(tx_tl_comp))

cat("\n— TX+TL + replication + septum (mother) —\n")
print(share(comp_with_fixed))

# Sentence for the manuscript:
mother_share <- (comp_with_fixed %>% dplyr::filter(compartment=="mother") %>% dplyr::pull(PT)) /
  sum(comp_with_fixed$PT)
foresp_share <- 1 - mother_share
cat(sprintf("\nSentence: Mother accounts for %.1f%% and the forespore for %.1f%% of total PT (TX+TL + replication + septum).\n",
            100*mother_share, 100*foresp_share))

#_______________________________________________________________________
# Share of energy used in the first hour (including genome replication)
# Manuscript: 
# Results: "Energetics of spore formation", second paragraph
#_______________________________________________________________________

# TX+TL spent in the first hour (0–1 h)
first_hour_T_TL <- hourly %>%
  filter(time_h == "1") %>%
  pull(total) %>%
  sum(na.rm = TRUE)

# Numerator: first-hour TX+TL + whole-genome replication
numerator   <- first_hour_T_TL + cost_rep_all
denominator <- all_pie_costs  # replication + TX + TL + septum

pct_first_hour <- 100 * numerator / denominator

# Compact report
summary_first_hour <- tibble::tibble(
  component = c("1st-hour TX+TL", "Genome replication", "Numerator total", "Total sporulation (denominator)"),
  PT        = c(first_hour_T_TL,    cost_rep_all,        numerator,         denominator),
  pretty    = sci(PT, 1)
)

print(summary_first_hour, n = Inf)
cat(sprintf(
  "\nApproximately %.1f%% of the total energy expenditure (including genome replication) occurs within the first hour of development.\n",
  round(pct_first_hour, 1)
))

#_____________________________________________________________________________
# How much energy is consumed until commitment (~2 h)?
# Assume molecules are recycled, how much of it non-refundable (direct costs)?
# Manuscript: 
# Results: "Energetics of spore formation", third paragraph 
#_____________________________________________________________________________
t_commit <- 2
hourly_num <- hourly %>% mutate(t = as.numeric(time_h))

# totals
tot_all <- sum(hourly_num$total, na.rm = TRUE)

# cumulative to commitment (≤2 h)
slice <- hourly_num %>% filter(t <= t_commit)

spent_to_commit <- sum(slice$total, na.rm = TRUE)
frac_spent_commit <- spent_to_commit / tot_all                     

direct_to_commit <- sum(slice$trscr_direct + slice$transl_direct, na.rm = TRUE)
frac_nonrecov_of_slice <- direct_to_commit / spent_to_commit       
frac_nonrecov_of_total <- direct_to_commit / tot_all               

cat(sprintf("By %0.1f h, incurred = %0.1f%% of total T+TL.\n",
            t_commit, 100*frac_spent_commit))
cat(sprintf("If precursors fully recycled: non-recoverable portion by %0.1f h = %0.1f%% of that slice (~%0.1f%% of total).\n",
            t_commit, 100*frac_nonrecov_of_slice, 100*frac_nonrecov_of_total))


                          ####################
########################## SPORE REVIVAL COSTS ##############################
                          ####################

#---------------------------------------------------------------------
# Replication costs of expressed genes during spore revival 
# Relative values compared to the whole genome
# Manuscript: 
# Results: "Energetics of spore revival", first paragraph
# Same equations as the whole genome synthesis, Eq.s 5-6
#----------------------------------------------------------------------

# Totals and relative fractions (length- and gene-count–based)
germ_related_genes_length <- sum(germ_once$gene_length, na.rm = TRUE)
germ_related_genes_number <- nrow(germ_once)
total_bacillus_genes      <- 4237  # same reference as in your spore block

rel_genome_germ <- germ_related_genes_length / genome.size
rel_genes_germ  <- germ_related_genes_number / total_bacillus_genes

# Replication costs for the germination/outgrowth gene set (same formula as sporeRep)
germRep <- germ_once %>%
  dplyr::mutate(
    opportunity = 2 * gene_length * dna_nucleotide_cost_opp,
    direct      = 2 * gene_length * dna_nucleotide_cost_dir +
      gene_length * helicase_cost +
      gene_length * primers_cost
  )

germRepSum <- dplyr::summarise(germRep,
                               sumOpp = sum(opportunity, na.rm = TRUE),
                               sumDir = sum(direct,      na.rm = TRUE))

germRepTotal <- with(germRepSum, sumOpp + sumDir)
germRepTotal_frac_of_genome <- germRepTotal / genome_tot # 11%

#_________________________________________________________________
# Revival Transcription
# Results: "Energetics of spore revival", first paragraph, 
# Figure 2 (PD, PO bars)
# Methods: "Bioenergetic accounting: definitions and assumptions, 
# "Transcription costs" 
# Equations: 7-9
#______________________________________________________________

# Transcription per interval (first; uses same logic as spores)
# - one-time synthesis + opportunity booked at first interval where w>0 (C_TD (synth) and C_TO)
# - polymerization (+ small supercoiling term) distributed by w (C_TD(t))

rna_nucleotide_cost_synthesis <- rna_nucleotide_cost_dir - poly_costs 

interval_order <- tibble(interval = interval_meta$interval,
                                 ord = seq_along(interval_meta$interval))

germ_tx <- germ_w %>%
  left_join(interval_order, by = "interval") %>%
  arrange(protID, ord) %>%
  group_by(protID) %>%
  mutate(
    first_hit        = which.max(w > 0), # first interval index with w>0
    row_in_group     = row_number(),
    is_first         = (row_in_group == first_hit),
    
    transcripts_total = proteins_per_cell / yield,
    transcripts_int   = transcripts_total * w,
    
    nt_total = gene_length * transcripts_total, # one-time
    nt_int   = gene_length * transcripts_int    # interval share
  ) %>%
  ungroup() %>%
  mutate(
    dir_synth_once = if_else(is_first, nt_total * rna_nucleotide_cost_synthesis, 0), # (c_sr)
    opp_once       = if_else(is_first, nt_total * rna_nucleotide_cost_opp, 0), # (c_or)
    
    dir_poly_int   = nt_int * poly_costs, # (c_pr)
    dir_super_int  = nt_int * c_supercoil # (c_sc)
  ) %>%
  group_by(interval, dur_h, mid_h) %>%
  summarise(
    transcription_dir = sum(dir_synth_once + dir_poly_int + dir_super_int, na.rm = TRUE),
    transcription_opp = sum(opp_once, na.rm = TRUE),
    .groups = "drop"
  )

#_______________________________________________________________
# Revival Translation
# Manuscript: 
# Results: "Energetics of spore revival", first paragraph, 
# Figure 2 (PD, PO bars)
# Methods: "Bioenergetic accounting: definitions and assumptions, 
# "Translational costs" 
# Equations: 10-11
#_______________________________________________________________

# Translation per interval (second) + modest extras (P_D only)
# (same extras toggles/params used for spores)
# extras can be switched on and off, minor contribution

# NA-safe AA cost fills
med_aa_dir <- median(germ_w$aa_directSum, na.rm = TRUE)
med_aa_opp <- median(germ_w$aa_opportunitySum, na.rm = TRUE)

# Per-interval translation (base + init/term only)
translation_interval <- germ_w %>%
  mutate(
    copies_int = proteins_per_cell * w,
    aa_dir_f   = coalesce(aa_directSum, med_aa_dir),
    aa_opp_f   = coalesce(aa_opportunitySum, med_aa_opp),
    
    # base translation costs
    direct_base = copies_int * aa_dir_f, # (C_Pr_D)
    opp_base    = copies_int * aa_opp_f, # (C_Pr_O)
    
    # extras (P_D): init + termination per protein copy
    pd_init_term = if (isTRUE(add_init_term)) copies_int * init_term_PD_per_protein else 0,
    
    # adjusted direct total
    direct_total = direct_base + pd_init_term
  ) %>%
  group_by(interval, dur_h, mid_h) %>%
  summarise(
    translation_dir = sum(direct_total, na.rm = TRUE),
    translation_opp = sum(opp_base, na.rm = TRUE),
    .groups = "drop"
  )

# Combine with transcription intervals and build totals
# First 15 min is germination #757214903 PT

germ_interval_costs <- translation_interval %>%
  full_join(germ_tx, by = c("interval","dur_h","mid_h")) %>%
  replace_na(
    list(translation_dir = 0, translation_opp = 0,
         transcription_dir = 0, transcription_opp = 0)
  ) %>%
  mutate(
    direct      = translation_dir + transcription_dir,
    opportunity = translation_opp + transcription_opp,
    total       = direct + opportunity
  ) %>%
  arrange(mid_h)

germ_totals <- germ_interval_costs %>%
  summarise(
    total_direct      = sum(direct),      # 1226210841 PD
    total_opportunity = sum(opportunity), # 4904910644 PO
    total_all         = sum(total)        # 6131121485 PT
  )

#___________________________________________________________
# Partial membrane synthesis
# Manuscript:
# Results: "Energetics of spore revival", first paragraph, 
# Figure 2 (pie chart)
# Methods: "Membrane synthesis and remodelling"
# Equations: 16
#___________________________________________________________

# knob
f_mem_rev <- 0.30   # fraction of the *whole-cell* membrane made during 0–3.5 h revival
f_recycle <- 1/6    # recycled fraction (i.e., 16.7% of lipid cost is avoided)

# safety clamp
f_mem_rev <- pmin(pmax(f_mem_rev, 0), 1)
f_recycle <- pmin(pmax(f_recycle, 0), 1)

# apply to whole-membrane totals
membraneRevOpp <- membraneOpp * f_mem_rev * (1 - f_recycle)
membraneRevDir <- membraneDir * f_mem_rev * (1 - f_recycle)
membraneRevTot <- membraneRevOpp + membraneRevDir 
# 54069922 PD, 636823528 PO, 690893450 PT 

#______________________________
# Total costs, plots, pie, bars
# Following will yield Figure 2
#______________________________

# Because the intervals are not the same stacked bars need some tweaks for aesthetics
# Bar sizes changed to reflect the reality but look nicer
# Backdrop totals (grey) and stacked PD/PO (colored)

# Pie components
pie_df <- germ_interval_costs %>%
  summarise(
    transcription = sum(transcription_dir + transcription_opp, na.rm = TRUE),
    translation   = sum(translation_dir   + translation_opp,   na.rm = TRUE)
  ) %>%
  pivot_longer(everything(), names_to = "category", values_to = "value")

# Include membraneGermTot
pie_df <- bind_rows(
  pie_df,
  tibble(category = "membrane", value = membraneRevTot) 
)

pie_df <- pie_df %>%
  mutate(share = 100 * value / sum(value),
         label = paste0(category, " ", sprintf("%.1f%%", share)))


pie_germ <- ggplot(pie_df, aes(x = "", y = value, fill = category)) +
  geom_col(width = 1) +
  coord_polar("y") +
  # mytheme +
  scale_fill_manual(values = c("transcription" = "#D55E00",
                               "translation"   = "#0072B2",
                               "membrane"      = "#CC79A7")) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  theme_void() +
  theme(legend.position = "none")

#ggsave("bioaccounting/figures/figure2_pieGerm.pdf", pie_germ, height = 5, width = 5)

## Time-series (stacked bars + asymptotic fit) 
bg_tot <- germ_interval_costs %>%
  transmute(time = mid_h, width = dur_h, total = total)

stack_df <- germ_interval_costs %>%
  select(mid_h, dur_h, opportunity, direct) %>%
  pivot_longer(c(opportunity, direct),
                      names_to = "type", values_to = "costs") %>%
  transmute(time = mid_h, width = dur_h, type, costs)

# Fit exp-decay to H0.5–H3.5 only (exclude the very first and the last bin) 
fit_df <- germ_interval_costs %>%
  filter(mid_h >= 0.375, mid_h <= 3.0) %>%   # 0.375=H0.5 mid, 3.0=H3.5 mid
  transmute(time = mid_h, costs = total)

fit_rev  <- nls(costs ~ SSasymp(time, yf, y0, log_alpha), data = fit_df)
pars     <- coef(fit_rev)
lambda_g <- exp(pars[["log_alpha"]])           # report as −λ in the label
pred_t   <- seq(min(fit_df$time), max(fit_df$time), length.out = 200)
pred_y   <- predict(fit_rev, list(time = pred_t))

# plain R^2 on linear scale (robust & readable)
RSSg <- sum((fit_df$costs - predict(fit_rev))^2)
TSSg <- sum((fit_df$costs - mean(fit_df$costs))^2)
R2_g <- 1 - RSSg/TSSg

pred_df <- data.frame(time = pred_t, y = pred_y)

# Plot
cols_PD_PO <- c(direct = "#D55E00", opportunity = "#0072B2")

pos <- position_dodge2(preserve = "single", padding = 0.05)

# precompute pretty labels using plotmath
lab_lambda <- sprintf("-lambda==%s", formatC(lambda_g, format = "f", digits = 3))
lab_r2     <- sprintf("R^2==%s",       formatC(R2_g,     format = "f", digits = 3))


f2 <- ggplot() +
  geom_col(data = bg_tot, aes(x = time, y = total, width = width),
           fill = "grey75", color = "grey90", alpha = 0.55) +
  geom_col(data = stack_df,
           aes(x = time, y = costs, fill = type, width = 0.8 * width),
           position = pos, color = "grey25") +
  geom_line(data = pred_df, aes(time, y), linewidth = 1) +
  scale_fill_manual(values = c(direct = "#D55E00", opportunity = "#0072B2"),
                    labels = c(direct = expression(P[D]),
                               opportunity = expression(P[O]))) +
  labs(x = "Time (h)", y = "ATP molecules") +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank())+
  annotate("text",
           x = max(pred_df$time), y = max(pred_df$y) * 0.20,
           label = lab_lambda, parse = TRUE,
           hjust = 1, size = 5) +
  annotate("text",
           x = max(pred_df$time), y = max(pred_df$y) * 0.14,
           label = lab_r2, parse = TRUE,
           hjust = 1, size = 5)+
 mytheme+
  theme(legend.position = c(0.3,0.9))+
  scale_y_continuous(breaks = c(4e8, 12e8, 20e8, 28e8), 
                     labels = c(4, 12, 20, 28), sec.axis=dup_axis())+
  scale_x_continuous(limits = c(0, 3.5), breaks = c(0, 1, 2, 3), labels = c(0, 1, 2, 3), sec.axis=dup_axis())+
  coord_cartesian(clip = "off") +
  annotate("text", x= 0.2, y = 3e9,label = paste("(x10^8)"), parse =T, size = 18/.pt)+
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 30))+
  annotate("text",
           x = 0.3, y = 12e8,
           label = "Germination", 
           vjust = -1, size = 5, fontface = "italic", angle = 90)

#ggsave("bioaccounting/figures/figure2_germCostsTime.pdf", f2, height = 5, width = 6)

#-------------------------------------------
# Practical summary for the manuscript text
#-------------------------------------------

ultimate_total_germ <- sum(pie_df$value) # 6823638259 PT with estimated membrane costs

# Slice intervals
germ_bins <- c("H0.25")
outg_bins <- c("H0.5","H1","H1.5","H2.5","H3.5")

sum_interval <- function(df, bins) {
  d <- df %>% dplyr::filter(interval %in% bins)
  tibble::tibble(
    PD   = sum(d$direct,            na.rm = TRUE),
    PO   = sum(d$opportunity,       na.rm = TRUE),
    PT   = PD + PO,
    TX_PT= sum(d$transcription_dir + d$transcription_opp, na.rm = TRUE),
    TL_PT= sum(d$translation_dir   + d$translation_opp,   na.rm = TRUE),
    TX_PD= sum(d$transcription_dir, na.rm = TRUE),
    TL_PD= sum(d$translation_dir,   na.rm = TRUE),
    TX_PO= sum(d$transcription_opp, na.rm = TRUE),
    TL_PO= sum(d$translation_opp,   na.rm = TRUE)
  )
}

germ_015  <- sum_interval(germ_interval_costs, germ_bins)
outgrowth <- sum_interval(germ_interval_costs, outg_bins)

# Revival membrane term (prefer *Rev*, else fallback to *Germ*)
mem_dir <- if (exists("membraneRevDir")) membraneRevDir else if (exists("membraneGermDir")) membraneGermDir else 0
mem_opp <- if (exists("membraneRevOpp")) membraneRevOpp else if (exists("membraneGermOpp")) membraneGermOpp else 0
mem_tot <- mem_dir + mem_opp

# Totals
revival_tx_tl_PT <- germ_015$PT + outgrowth$PT
revival_total_PT <- revival_tx_tl_PT + mem_tot

# Shares within REVIVAL
share_tx <- 100 * (germ_015$TX_PT + outgrowth$TX_PT) / revival_total_PT
share_tl <- 100 * (germ_015$TL_PT + outgrowth$TL_PT) / revival_total_PT
share_mm <- 100 *  mem_tot / revival_total_PT

# Spore formation totals you already have
spore_T_PT  <- (if (exists("cost_transcript")) cost_transcript else sum(transcription_hourly$total, na.rm = TRUE)) +
  (if (exists("cost_translation")) cost_translation else sum(translation_by_time_extras$total, na.rm = TRUE))
spore_T_plus_septum_PT <- spore_T_PT + (if (exists("septumTot")) septumTot else 0)

# Pretty print tables
ms_rev_summary <- tibble::tribble(
  ~item,                                   ~PT,                    ~PD,             ~PO,            ~note,
  "Germination (0–0.25 h)",                germ_015$PT,           germ_015$PD,     germ_015$PO,    "TX+TL only",
  "Outgrowth (0.5–3.5 h)",                 outgrowth$PT,          outgrowth$PD,    outgrowth$PO,   "TX+TL only",
  "Membrane during revival",               mem_tot,               mem_dir,         mem_opp,        "with recycling factor",
  "Revival total (TX+TL+membrane)",        revival_total_PT,      NA,              NA,             "sum of rows above",
) %>%
  dplyr::mutate(
    PT_str = sci(PT, 1),
    PD_str = dplyr::if_else(is.na(PD), "", sci(PD, 1)),
    PO_str = dplyr::if_else(is.na(PO), "", sci(PO, 1))
  ) %>%
  dplyr::select(item, PT_str, PD_str, PO_str, note)

revival_shares <- tibble::tibble(
  category = c("transcription", "translation", "membrane"),
  percent  = c(share_tx, share_tl, share_mm)
)

# Helper manuscript text 
cat("\n=== Revival summary (PT, PD, PO) ===\n")
print(ms_rev_summary, n = Inf)
cat("\n=== Shares within revival total (%) ===\n")
print(revival_shares)

cat("\nProse helpers:\n")
cat(sprintf("\nUltimate revival total (TX + TL + membrane): %s PT\n",
            sci(ultimate_total_germ, 1)))
cat(sprintf("Revival total (TX+TL+membrane): %s PT\n", sci(revival_total_PT)))
cat(sprintf("Germination 0–0.25 h: %s PT (%s PD, %s PO)\n",
            sci(germ_015$PT), sci(germ_015$PD), sci(germ_015$PO)))
cat(sprintf("Outgrowth 0.5–3.5 h: %s PT (%s PD, %s PO)\n",
            sci(outgrowth$PT), sci(outgrowth$PD), sci(outgrowth$PO)))
cat(sprintf("Membrane during revival: %s PT (%s PD, %s PO)\n",
            sci(mem_tot), sci(mem_dir), sci(mem_opp)))
cat(sprintf("Revival composition: TX %.1f%% | TL %.1f%% | Membrane %.1f%%\n",
            share_tx, share_tl, share_mm))

####################################
# Head to Head comparisons of traits 
# Following will yield Figure 3
# Manuscript: 
# Results: "Head-to-head comparison"
####################################

# Other cellular structures (non-sporulation/germination/outgrowth)
# Build-only costs per cell (no maintenance)

# constants (already defined upstream) 
# poly_costs <- 2
# rna_nucleotide_cost_synthesis <- rna_nucleotide_cost_dir - poly_costs  # c_sr
# yield <- 100  # proteins per mRNA lifetime

# Medians for safe fills 
# prot_len_med 
# gene_len_med

#  Map trait genes → locus_tag via SubtiWiki symbols (primary), then join protein table
#  Exclude developmental program categories to avoid double-counting with time-series blocks
trait_floor_copies <- 0  # set to >0 to increase floor

traits_base <- otherTraits %>%
  filter(!category %in% c("sporulation", "germination", "outgrowth")) %>%
  mutate(gene_std = std(gene)) %>%
  left_join(
    annotation_clean %>% select(gene_std, locus_tag, gene_length, protein_length),
    by = "gene_std"
  ) %>%
  # bring in abundance-derived fields and AA costs
  left_join(
    protSeqTidyAbun %>%
      select(locus_tag, proteins_per_cell, aa_directSum, aa_opportunitySum),
    by = "locus_tag") %>%
  # prefer numeric lengths from prot table; fall back to annotation; then to medians
  mutate(
    gene_length      = coalesce(gene_length,  gene_len_med),
    protein_length   = coalesce(protein_length, prot_len_med),
    proteins_per_cell = coalesce(proteins_per_cell, trait_floor_copies)  # no imputation here; 0 if unseen
  ) %>%
  select(category, gene = gene, locus_tag, proteins_per_cell,
                gene_length, protein_length, aa_directSum, aa_opportunitySum)

# Fill missing AA cost sums conservatively with medians (per-protein totals)
med_aa_dir <- median(traits_base$aa_directSum,      na.rm = TRUE)
med_aa_opp <- median(traits_base$aa_opportunitySum, na.rm = TRUE)

traits_costed <- traits_base %>%
  mutate(
    aa_dir_f = coalesce(aa_directSum,      med_aa_dir),
    aa_opp_f = coalesce(aa_opportunitySum, med_aa_opp),
    
    # Translation build: copies × per-protein AA cost
    trans_direct = proteins_per_cell * aa_dir_f,
    trans_opp    = proteins_per_cell * aa_opp_f,
    trans_total  = trans_direct + trans_opp,
    
    # Transcription build: one-time synthesis of the mRNAs needed
    transcripts_total = proteins_per_cell / yield,             # copies per cell
    nt_total          = transcripts_total * gene_length,       # nucleotides to make those mRNAs
    tx_direct_once    = nt_total * rna_nucleotide_cost_synthesis,
    tx_opp_once       = nt_total * rna_nucleotide_cost_opp,
    tx_total_once     = tx_direct_once + tx_opp_once
  )

#  Aggregate per category (build-only; no per-gene replication here)
totalCosts_traits <- traits_costed %>%
 group_by(category) %>%
    summarise(
    translation_direct        = sum(trans_direct,   na.rm = TRUE),
    translation_opportunity   = sum(trans_opp,      na.rm = TRUE),
    transcription_direct      = sum(tx_direct_once, na.rm = TRUE),
    transcription_opportunity = sum(tx_opp_once,  na.rm = TRUE),
    total                     = translation_direct + translation_opportunity +
      transcription_direct + transcription_opportunity,
    .groups = "drop"
  ) %>%
  select(category, total)
  
# Build a top-level comparison table with global items
#    - Program totals computed above:
#        spore_total_all  = sum(transcription_hourly$total) + sum(translation_by_time_extras$total)
#        germ_total_all   = sum(germ_interval_costs$total) (first 15 min is germination)
#    - Use membrane cost from the membrane block (membraneTot)
#    - Use genome replication total (genome_tot)
#    - Use total cell budget C_T for a reference bar (or separate build vs maintenance)

## Constants (Lynch & Marinov 2015, pulled from the supplementary table; 20 °C standard)
C_M_per_h  <- 1.159e9    # ATP / cell / hour (maintenance)
C_G_cell   <- 9.251e10  # ATP / cell (growth/build)
gen_time_h <- 1.16      # h generation at 20 °C 
budget_ref <- C_G_cell + gen_time_h * C_M_per_h  # "total cell budget" bar

## Program window for maintenance 
# Set sporulation hours from your dataset; if unknown, keep a sensible default
t_spor <- 8 # hours
t_germ <- 0.25 # first 15 min
t_outg <- 3.25 # Until 210 min 
t_program <- t_spor + t_germ + t_outg

maintenance_prog <- C_M_per_h * t_program  # ATP

## Developmental program total (build-only) 
cost_sporulation <- cost_transcript + cost_translation

cost_germ_015 <- germ_interval_costs %>%
  filter(interval == "H0.25") %>%
  summarise(v = sum(total, na.rm = TRUE)) %>% pull(v)

cost_outgrowth_rest <- germ_interval_costs %>%
  filter(interval %in% c("H0.5","H1","H1.5","H2.5","H3.5")) %>%
  summarise(v = sum(total, na.rm = TRUE)) %>% pull(v)

dev_parts <- tibble(
  phase = c("germination (0.25 h)", "spore formation", "outgrowth"),
  total = c(cost_germ_015, cost_sporulation,  cost_outgrowth_rest)
)

dev_total <- sum(cost_sporulation, cost_germ_015, cost_outgrowth_rest) 

## Trait totals you already computed (build-only) 
# totalCosts_traits : data frame with columns category, total  (build-only by category)
# Add your whole-membrane & genome replication totals if available

## Build the bars table (all grey, no stacks)
bars_df <- totalCosts_traits %>%
  select(category, total) %>%
  add_row(category = "membrane lipid synthesis",
                 total    = if (exists("membraneTot")) membraneTot else NA_real_) %>%
  add_row(category = "genome replication",
                 total    = if (exists("genome_tot"))  genome_tot  else NA_real_) %>%
  add_row(category = "developmental program", total = dev_total) %>%
  add_row(category = "maintenance (program window)", total = maintenance_prog) %>%
  add_row(category = "total cell budget", total = budget_ref) %>%  # anchor for %
  filter(is.finite(total), total > 0)

pretty_names <- c(
  "membrane lipid synthesis"     = "Membrane lipids",
  "genome replication"           = "Genome replication",
  "developmental program"        = "Spore life cycle",
  "maintenance (program window)" = "Maintenance during cycle",
  "total cell budget"            = "Total cell budget",
  "swarming"                     = "Swarming", 
  "chemotaxis_motility"          = "Chemotaxis", 
  "biofilm"                      = "Biofilm", 
  "flagella"                     = "Flagella", 
  "essential"                    = "Essential genes", 
  "competence"                   = "Competence", 
  "heat_shock"                   = "Heat shock proteins", 
  "homeostasis"                  = "Homeostasis"
  
)

# keep only keys that actually exist to avoid warnings
pretty_names <- pretty_names[names(pretty_names) %in% unique(bars_df$category)]

bars_nice <- bars_df %>%
  arrange(total) %>%                                
  mutate(
    category_pretty = recode(category, !!!pretty_names, .default = category),
    category_pretty = factor(category_pretty, levels = unique(category_pretty))
  )

# Cumulative dotted marks for the program bar 
dev_name_pretty <- recode("developmental program", !!!pretty_names)

y_dev <- match(dev_name_pretty, levels(bars_nice$category_pretty))

seg_df <- tibble(
  x    = cumsum(dev_parts$total),  # germ, germ+spore, germ+spore+outgrowth
  y    = y_dev,
  ymin = y_dev - 0.38,
  ymax = y_dev + 0.38
)

pct_breaks <- budget_ref * c(1, 10, 100) / 100  # 1%, 10%, 100% of budget

xmin <- min(bars_nice$total[bars_nice$total > 0]) * 0.8
xmax <- max(bars_nice$total) * 1.25

# Fancy scientific axis names
scientific_10 <- function(x) {
  raw <- ifelse(x == 0, "0",
                sub("e[+]?", " %*% 10^", scales::scientific_format()(x)))
  # drop "1 %*% 10^" → "10^"
  neat <- gsub("^1\\s*%\\*%\\s*10\\^", "10^", raw)
  parse(text = neat)
}

bar_width <- 0.8  # was 0.8 → more vertical gap between categories
tick_half <- bar_width/2 + 0.1

f3 <- ggplot() +
  geom_col(
    data = bars_nice,
    aes(x = total, y = category_pretty),
    width = bar_width, fill = "grey80", color = "black"
  ) +
  scale_x_log10(
    name   = "ATP molecules",
    labels = scientific_10,
    expand = expansion(mult = c(0.08, 0.03)),            # left gutter
    sec.axis = sec_axis(~ ., breaks = pct_breaks, labels = c("1","10","100"),
                        name = "Costs relative to the cell budget (%)")
  ) +
  coord_cartesian(xlim = c(xmin * 0.8, xmax)) +          # keep default clip = "on"
  mytheme +
  theme(
    axis.title.x.top   = element_text(hjust = 0.5, size = 18),
    axis.text.x.top    = element_text(size = 16),
    axis.text.y        = element_text(margin = margin(r = 8)),
    panel.grid.major.y = element_blank(),                # prevent “strikethrough”
    panel.ontop        = FALSE                           # just in case mytheme sets it TRUE
  ) +
  ylab(NULL) +
  geom_segment(
    data = seg_df |> dplyr::mutate(ymin = y - tick_half, ymax = y + tick_half),
    aes(x = x, xend = x, y = ymin, yend = ymax),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7
  )

#ggsave("bioaccounting/figures/figure3_comparison.pdf", f3, height = 6, width = 7.5)


#__________________________
# Numbers for text/caption 
#__________________________

# Table with % of total cell budget
head_to_head_tbl <- bars_nice %>%
  transmute(
    category = as.character(category_pretty),
    ATP      = total,
    pct_of_budget = 100 * total / budget_ref
  ) %>%
  arrange(desc(ATP))

print(head_to_head_tbl, n = Inf)

# Verify the dotted ticks equal the cumulative phase totals
dev_bar_total <- head_to_head_tbl$ATP[head_to_head_tbl$category == dev_name_pretty]

# One-liners
cat(sprintf("\nSpore life cycle = %.2f%% of total cell budget\n",
            100 * dev_bar_total / budget_ref))
cat(sprintf("Maintenance during cycle = %.2f%% of total cell budget\n",
            100 * head_to_head_tbl$ATP[head_to_head_tbl$category == 'Maintenance during cycle'] / budget_ref))
cat(sprintf("Genome replication = %.2f%% | Membrane lipids = %.2f%% of budget\n",
            100 * head_to_head_tbl$ATP[head_to_head_tbl$category == 'Genome replication'] / budget_ref,
            100 * head_to_head_tbl$ATP[head_to_head_tbl$category == 'Membrane lipids'] / budget_ref))

#________________________________________________________________________________
# Spore life cycle energetic investment relative to total cellular energy budget.
# Here all costs are included: build costs + genome + membrane
# Manuscript: 
# Results, first (intro) paragraph
#________________________________________________________________________________

all_cycle <- (dev_total + genome_tot + membraneRevTot + septumTot)
relative_cost <- all_cycle/budget_ref #10% 

#________________________________________________________________________________
# REVISION - 19 November 2025 
#________________________________________________________________________________

####################################
# Figure 3: Per-generation-time normalized costs
# Revised to address reviewer comments
####################################

# Define durations for traits that span multiple generations
# NA means "instantaneous" or "per-generation already"
trait_durations <- tribble(
  ~category,              ~duration_h,  ~notes,
  "biofilm",              12,           "Time to robust pellicle formation",
  "competence",           3,            "Transient K-state duration",
  "flagella",             NA_real_,     "Persistent structure, expressed per generation",
  "swarming",             NA_real_,     "Uses flagella; expressed per generation",
  "chemotaxis_motility",  NA_real_,     "Instantaneous - already per generation",
  "essential",            NA_real_,     "Instantaneous - already per generation", 
  "heat_shock",           NA_real_,     "Instantaneous - already per generation",
  "homeostasis",          NA_real_,     "Instantaneous - already per generation"
)

## Constants (Lynch & Marinov 2015)
C_M_per_h  <- 1.159e9    # ATP / cell / hour (basal maintenance)
C_G_cell   <- 9.251e10   # ATP / cell (growth/build)
gen_time_h <- 1.16       # generation time at 20°C (hours)
budget_ref <- C_G_cell + gen_time_h * C_M_per_h  # Total cell budget per generation

## Normalize all costs to per vegetative generation (1.16 h)
# For programs: divide by (duration_h / gen_time_h)
# For instantaneous traits: no adjustment (duration_h = NA or gen_time_h)

totalCosts_normalized <- totalCosts_traits %>%
  left_join(trait_durations, by = "category") %>%
  mutate(
    # Calculate how many generation-times this trait spans
    generation_times = coalesce(duration_h / gen_time_h, 1.0),
    
    # Normalize to per generation
    # If instantaneous (NA duration), generation_times = 1, no adjustment
    # If multi-generation program, divide by generation_times
    total = total/ generation_times
  ) %>%
  select(category, total, duration_h, generation_times)

# Report normalization
cat("\n=== Trait normalization to per vegetative generation ===\n")
totalCosts_normalized %>%
  transmute(
    Category = category,
    Duration_h = coalesce(duration_h, gen_time_h),
    Gen_times = round(generation_times, 2),
    Normalized_cost = scales::scientific(total, digits = 2),
    Adjustment = ifelse(generation_times > 1, 
                        sprintf("÷%.1f", generation_times), 
                        "none")
  ) %>%
  print(n = Inf)

## Developmental program: ONE complete spore generation (11.5 h)
# This is NOT normalized - it represents one biological generation
# for the spore life history pathway
t_spor <- 8      # spore formation (hours)
t_germ <- 0.25   # germination (first 15 min)
t_outg <- 3.25   # outgrowth (until 210 min)
t_program <- t_spor + t_germ + t_outg  # 11.5 h total

# Spore cycle costs (transcription + translation only for comparison)
cost_sporulation <- cost_transcript + cost_translation

cost_germ_015 <- germ_interval_costs %>%
  filter(interval == "H0.25") %>%
  summarise(v = sum(total, na.rm = TRUE)) %>% pull(v)

cost_outgrowth_rest <- germ_interval_costs %>%
  filter(interval %in% c("H0.5","H1","H1.5","H2.5","H3.5")) %>%
  summarise(v = sum(total, na.rm = TRUE)) %>% pull(v)

dev_total <- sum(cost_sporulation, cost_germ_015, cost_outgrowth_rest)

## Maintenance per generation-time (NOT per program window)
maintenance_per_gentime <- C_M_per_h * gen_time_h

## Build the bars table
# Genome & membrane: already per vegetative generation (replication during 1.16h cycle)
bars_df <- totalCosts_normalized %>%
  select(category, total) %>%
  add_row(category = "membrane lipid synthesis", 
          total = membraneTot) %>%  # whole cell membrane, per veg gen
  add_row(category = "genome replication",       
          total = genome_tot) %>%   # whole genome, per veg gen
  add_row(category = "developmental program",    
          total = dev_total) %>%    # ONE complete spore generation (11.5h)
  add_row(category = "maintenance",              
          total = maintenance_per_gentime) %>%  # per veg gen (1.16h)
  add_row(category = "total cell budget",        
          total = budget_ref) %>%   # reference (C_G + C_M*gen_time)
  filter(is.finite(total), total > 0)

# Pretty names for display
pretty_names <- c(
  "membrane lipid synthesis"     = "Membrane lipids",
  "genome replication"           = "Genome replication",
  "developmental program"        = "Spore life cycle",
  "maintenance"                  = "Maintenance energy",
  "total cell budget"            = "Total cell budget",
  "swarming"                     = "Swarming", 
  "chemotaxis_motility"          = "Chemotaxis", 
  "biofilm"                      = "Biofilm", 
  "flagella"                     = "Flagellum", 
  "essential"                    = "Essential genes", 
  "competence"                   = "Competence", 
  "heat_shock"                   = "Heat shock", 
  "homeostasis"                  = "Homeostasis"
)

bars_nice <- bars_df %>%
  arrange(total) %>%
  mutate(
    category_pretty = recode(category, !!!pretty_names, .default = category),
    category_pretty = factor(category_pretty, levels = unique(category_pretty))
  )

# For the spore bar, add cumulative phase markers
dev_name_pretty <- recode("developmental program", !!!pretty_names)
y_dev <- match(dev_name_pretty, levels(bars_nice$category_pretty))

dev_parts <- tibble(
  phase = c("germination (0.25 h)", "spore formation", "outgrowth"),
  total = c(cost_germ_015, cost_sporulation, cost_outgrowth_rest)
)

seg_df <- tibble(
  x    = cumsum(dev_parts$total),
  y    = y_dev,
  ymin = y_dev - 0.38,
  ymax = y_dev + 0.38
)

# Plot setup
pct_breaks <- budget_ref * c(0.1, 1, 10, 100) / 100
xmin <- min(bars_nice$total[bars_nice$total > 0]) * 0.8
xmax <- max(bars_nice$total) * 1.25
primary_breaks <- 10^(8:11)
primary_breaks <- primary_breaks[primary_breaks >= xmin & primary_breaks <= xmax]

pct_labeller <- function(x) {
  sapply(x, function(val) {
    diffs <- abs(log10(val) - log10(pct_breaks))
    if(min(diffs) < 0.05) {
      idx <- which.min(diffs)
      return(c("0.1", "1", "10", "100")[idx])
    } else {
      return("")
    }
  })
}

bar_width <- 0.8
tick_half <- bar_width/2 + 0.1

# Create plot
f3 <- ggplot() +
  geom_col(
    data = bars_nice,
    aes(x = total, y = category_pretty),
    width = bar_width, fill = "grey80", color = "black"
  ) +
  scale_x_log10(
    name   = "ATP molecules per generation",
    breaks = primary_breaks,
    labels = trans_format("log10", math_format(10^.x)),
    expand = expansion(mult = c(0.08, 0.03)),
    sec.axis = sec_axis(~ ., 
                        breaks = pct_breaks, 
                        labels = pct_labeller,
                        name = "Costs relative to the cell budget (%)")
  ) +
  coord_cartesian(xlim = c(xmin * 0.8, xmax)) +
  mytheme +
  theme(
    axis.title.x.top   = element_text(hjust = 0.5, size = 18),
    axis.text.x.top    = element_text(size = 16),
    axis.text.y        = element_text(margin = margin(r = 8)),
    panel.grid.major.y = element_blank(),
    panel.ontop        = FALSE
  ) +
  ylab(NULL) +
  geom_segment(
    data = seg_df %>% mutate(ymin = y - tick_half, ymax = y + tick_half),
    aes(x = x, xend = x, y = ymin, yend = ymax),
    inherit.aes = FALSE, linetype = "dotted", linewidth = 0.7
  )

ggsave("bioaccounting/figures/figure3_comparison_revised.pdf", f3, height = 6, width = 7)

# Summary table for manuscript
head_to_head_tbl <- bars_nice %>%
  transmute(
    category = as.character(category_pretty),
    ATP      = total,
    pct_of_budget = 100 * total / budget_ref
  ) %>%
  arrange(desc(ATP))

print(head_to_head_tbl, n = Inf)

cat(sprintf("\n=== Normalization summary ===\n"))
cat(sprintf("Spore life cycle: %.1f h (%.1f veg gen-times) = 1 spore generation\n",
            t_program, t_program/gen_time_h))
cat(sprintf("All other traits: Normalized to per vegetative generation (%.2f h)\n", gen_time_h))
cat(sprintf("  - Biofilm: %.0f h ÷ %.1f gen = per-gen cost\n", 12, 12/gen_time_h))
cat(sprintf("  - Competence: %.0f h ÷ %.1f gen = per-gen cost\n", 3, 3/gen_time_h))
cat(sprintf("  - Instantaneous traits: No adjustment\n\n"))