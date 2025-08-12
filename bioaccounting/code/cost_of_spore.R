# Spore Costs
# 21 July 2025 - last update
# Author: C. Karakoc
# Global expression data from SporeWeb: https://sporeweb.molgenrug.nl/
# Newly synthesized proteins during germination: doi: https://10.1128/mSphere.00463-20
# Global protein abundance data: https://pax-db.org/
# List of gene categories & annotation: https://subtiwiki.uni-goettingen.de/
# Protein sequence: Uniprot https://www.uniprot.org/taxonomy/224308
# Amino acid costs: https://doi.org/10.1073/pnas.1701670114
# COGs: https://journals.asm.org/doi/10.1128/jb.00079-22

######################
# Packages & Plotting 
######################
##########################################################################
library(tidyverse)       
library(patchwork) # for nls
library(stringr) 
#devtools::install_github("dgrtwo/fuzzyjoin")
library(fuzzyjoin) # for merging protein sequences
# for COGs
library(vegan)
# ggplot theme
library(scales) 
library(ggfortify)
library(ggrepel)
#plotting
library(ggpubr) 
library(grid)

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

mytheme2 <- theme_bw()+
  theme(axis.ticks.length = unit(.25, "cm"))+
  theme(legend.text = element_text(size=14))+
  theme(axis.text = element_text(size = 18, color = "black"), axis.title = element_text(size = 18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.border = element_rect(fill=NA, colour = "black", 
                                    size=1))+
  theme(strip.text.x = element_text(size = 18))+
  theme(legend.title=element_blank())+
  theme(panel.border = element_rect(fill=NA, colour = "black", 
                                    linewidth=1)) +
  #theme(axis.text.x.top = element_blank(), axis.title.x.top = element_blank(),
  #     axis.text.y.right = element_blank(), axis.title.y.right = element_blank())+
  theme(axis.title.x = element_text(margin=margin(10,0,0)),
        axis.title.y = element_text(margin=margin(0,10,0,0)),
        axis.text.x = element_text(margin=margin(10,0,0,0)),
        axis.text.y = element_text(margin=margin(0,10,0,0)))+
  theme(axis.title.x.top = element_text(margin=margin(10,10,0)),
        axis.title.y.right = element_text(margin=margin(0,10,0,0)),
        axis.text.x.top = element_text(margin=margin(10,0,10,0)),
        axis.text.y.right = element_text(margin=margin(0,10,0,0)))

# Color blind palette
cbpalette <- c("#0072B2", "#D55E00","#009E73", "#CC79A7", "#56B4E9", "#999999", "#F0E442", "#000000")

#setwd("~/Documents/GitHub/cost_of_spore")
setwd("~/GitHub/cost_of_spore")

######################
# Data
######################
##########################################################################
# Gene & protein length from SubtiWiki
annotationData <- read.table("./bioaccounting/data/subtiwiki.gene.export.2022-05-11.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Data with synonyms # different sources
nameMap        <- read.table("./bioaccounting/data/nameMap.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Protein abundances from PaxDB
protAbun       <- read.table("./bioaccounting/data/protAbunData.csv", sep = ',', dec = ".", header = T, stringsAsFactors = F, na.strings=c(" ","NA"))

# Protein sequences from Uniprot
protSeq        <- read.delim("./bioaccounting/data/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.02.08-14.23.07.59.txt")

# Nucleotide & Amino acid costs 
aaCosts        <- read.table("./bioaccounting/data/aaCosts_pnas.1701670114.sd01.csv", sep = ",", dec = "." , header = T)
nucCosts       <- read.table("./bioaccounting/data/nucleotideCosts_pnas.1701670114.sd01.csv", sep = ",", dec = "." , header = T)

# Other traits from SubtiWiki lists 
otherTraits    <- read.table("./bioaccounting_data/traits.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

# Germination time course - Swarge et al. 2020 
germination6   <- read.table("./bioaccounting/data/SwargeEtAl_onlyNewProteins.csv", sep = ",", header = T) 
  
# COGs - Galperin
meta_data   <- read.table("./bioaccounting/data/calc_meta.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

                            #####################
########################### SPORE FORMATION COSTS #####################################
                            #####################
#####################
# Time series prep
#####################

# Expression data from SporeWeb
files  <- list.files(path = "./bioaccounting/data/SporeWebHeatmaps/" , pattern = "*.csv", full.names = T)
exp_files <-  list()
for (i in 1:length(files)){
  exp_files[[i]] <- read.table(files[i], header = T, sep = ",", dec = ".")
}
merged_exp                   <- bind_rows(exp_files, .id = 'sourceID')
merged_exp$sourceID          <- as.factor(merged_exp$sourceID)
levels(merged_exp$sourceID ) <- c("1.vegetative", "2.starvation", "3.onset", "4.commitment", "5.engulfment")

# Wrangle 
expressionLong <- merged_exp[,c(1, 3:5, 7:14)] %>%
  pivot_longer(cols = c('t1','t2','t3','t4','t5','t6','t7','t8'),
               names_to =  "time", values_to = "expression") %>%
  group_by(time, sourceID, regulators, gene, locus_tag) %>%
  summarise(meanexp = mean(expression)) %>%
  ungroup() %>%
  filter(meanexp > 0) %>%
mutate(time_h = gsub('t', '', time)) 

# Abundance data prep
#####################
# Protein abundance
#####################
##########################################################################

# Protein abundance data does not include locus tags
# It is often easier to merge data sets with locus tags, because genes have a lot of synonyms

# A gene name assigned two different locus tag, corrected based on PaxDB
which(nameMap$gene1 =="nrgB")
nameMap[246,4] <- ""

# Logical map with matching gene names
MAP <- as.data.frame(t(apply(nameMap, 1, function(x) x %in% protAbun$gene)))
colnames(MAP) <- colnames(nameMap)

# Convert logical map into a column mapped with the matching column names
MAP$index = apply(MAP, 1, function(x) paste(names(x)[x], collapse=", "))
nameMap$index = MAP$index

# Match the column names by the rows
id <- lapply(seq(nrow(nameMap)), function(i) nameMap[i,nameMap$index[i]])

# Convert the list output to vector column, but prevent loosing NULLs, convert them
# First to NAs, so that we know which genes do not have abundance info
nameMap$gene <- do.call(rbind, lapply(id, function(x) if(is.null(x)) NA else x))
nameMap$gene <- ifelse(is.na(nameMap$gene), nameMap$geneP, nameMap$gene)

# Now use this collapsed matching column to merge abundances, genes and locus tag
# Merge data
mergedAbunData <- protAbun  %>%
  right_join(nameMap, by = "gene", multiple = "all") %>%
  left_join(annotationData[-2], by = "locus_tag", multiple = "all") %>%
  distinct(locus_tag, .keep_all = TRUE) %>%  #left join duplicates 
  dplyr::select(locus_tag, gene, abundance, protein_length, gene_length)

# Fill NAs with median values

# Median protein abundance
# Manuscript line: 1,045
protMed <- as.numeric(as.vector(protAbun$abundance))
median.protMed <- median(protMed, na.rm = T) #17.6 #round
median.protMed

mergedAbunData$protein_length <- as.numeric(mergedAbunData$protein_length)
mergedAbunData$gene_length    <- as.numeric(mergedAbunData$gene_length)

# Median protein and gene length 
median.protein_length <- median(mergedAbunData$protein_length, na.rm = T) #254
median.gene_length <- median(mergedAbunData$gene_length, na.rm = T) #765
median.protein_length
median.gene_length

# This function takes a while to run 
protSeqTidy <- protSeq %>%
  regex_inner_join(mergedAbunData, by = "gene") %>%
  distinct(sequence, .keep_all = TRUE) #remove the duplicated rows bases on unique sequences
 
# Amino acid alphabet 
alphabet = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

seqCount <- protSeqTidy %>%
  rowwise() %>%
  reframe(aac = str_count(sequence, pattern = alphabet)) %>%
  mutate(symbol = rep(alphabet, times = 4243)) %>%
  mutate(gene = rep(1:4243, each = 20)) %>%
  left_join(aaCosts, by = "symbol") %>%
  mutate(aa_opportunity = aac*opportunity_costs, aa_direct = aac*direct_costs) %>%
  group_by(gene) %>%
  summarize(aa_opportunitySum = sum(aa_opportunity), aa_directSum = sum(aa_direct)) %>%
  mutate(locus_tag = protSeqTidy$locus_tag) %>%
  distinct(locus_tag, .keep_all = TRUE)
  
protSeqTidyAbun <-  protSeqTidy %>%
  left_join(seqCount[,-1], by = "locus_tag") %>%
  select("gene.y", "protID", "locus_tag", "abundance", "protein_length", "gene_length", "aa_opportunitySum", "aa_directSum") %>%
  distinct(locus_tag, .keep_all = TRUE)


# Accounting data prep
#####################
# Expression data
#####################

# Merge with expression data
# Here I'll create two different data sets. One for calculating transcription
# costs, another for translation. Since proteins are degraded much slower, I
# I will account for only repolymerization costs of transcripts 
# Here genes are accounted once base on first appearance

mergedExpData_time_distinct <- expressionLong %>%
  distinct(locus_tag, .keep_all = TRUE) %>% # genes are accounted only once
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, median.protMed)) %>%
  mutate(gene_length.filled = replace_na(gene_length, median.gene_length)) %>%
  mutate(protein_length.filled = replace_na(protein_length, median.protein_length))

mergedExpData_time <- expressionLong %>%
  left_join(protSeqTidyAbun, by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, median.protMed)) %>%
  mutate(gene_length.filled = replace_na(gene_length, median.gene_length)) %>%
  mutate(protein_length.filled = replace_na(protein_length, median.protein_length))

# Replication costs (Whole genome) 
#####################
# Genome costs
#####################
##########################################################################
# Genome size = (https://www.nature.com/articles/36786)

# DNA unwinding 1 ATP per base pair #https://doi.org/10.1016/j.cell.2006.10.049
# ligation of Okazaki fragments 0.3ATP #Lynch and Marinov, 2015
# included 

# Opportunity 
# 2 * Genome size * (34+1)
# Eq. 6
genome.size <- 4214810
genome_opp <- 2 * genome.size * 35 #295036700
# Direct 
# 2 * Genome size * (11+2)
# Eq. 5
genome_dir <-( 2 * genome.size * 14 ) # 118014680
# Total
genome_tot <- genome_opp + genome_dir #413051380

#####################
# Membrane costs
#####################

# Number of lipid molecules = Cellular membrane areas/head-group areas of membrane lipid molecules
# Head group area is a1 = 0.65 nm2 (Nagle and Tristram-Nagle 2000; Petrache et al. 2000; Kucerka et al. 2011).

# Thickness of the bilayer (h):
# The thickness of a bilayer is approximately twice the radius of the head-group area, which 0.5 nm in all cases,
# plus the total length of the internal hydrophobic tail domains (Lewis and Engelman 1983; Mitra et al. 2004), 
# generally 3.0 nm, so total is 4 nm

# Bacillus average length (a) and width (b) Barak et al. 2018 
a1 <- 0.65
h  <- 4
a  <- 2.5*1000 # convert to nm
b  <- 1*1000 # also septum
# height / width = 2.5
c  <- 0.4 

# Outer # 4πa*b
# Eq. 11
outer <- 4*pi*a*b 
moleculesOut <- outer/a1 
# Inner membrane 4π(a − h)(b - h)
inner <- 4*pi*((a-h)*(b-h))
moleculesInn <- inner/a1 
# 50% discount, because of proteins
totalMol <- (moleculesOut + moleculesInn)/2 # total lipid molecules

# Cost of lipid head & tail 
# Opportunity = 212, Direct = 18
# costs
lipid.opp <- 212
lidid.dir <- 18
membraneOpp <- totalMol*lipid.opp 
membraneDir <- totalMol*lidid.dir 
membraneTot <- membraneOpp + membraneDir

# Septum should be 1µm 1000 nm (as the width of the cell) 
septumOut <- ((4*b)/a1)/2
septumInn <- ((4*(b-h))/a1)/2 
septumTot <- septumOut + septumInn

# costs
septumOpp     <- septumTot*lipid.opp 
septumDir     <- septumTot*lidid.dir 

# Germination assuming that they recycle membrane of the endospore
# Whole membrane - (endospore sphere + septum)
# Septum stretches 
# Endospore size is 1/6 of the total cell

membraneGerm     <- membraneTot - (membraneTot/6) #9237762011
membraneGermDir  <- membraneDir - (membraneDir/6) #722955288
membraneGermOpp  <- membraneOpp - (membraneOpp/6) #8514806723

# Replication costs of expressed genes
#####################
# Spore Replication 
#####################

# 878 genes
# 4429 genes whole genome
sum(mergedExpData_time_distinct$gene_length.filled)
# total length = 722878
# genome length = 4214810
# 722878/4214810 %17
# 878/4429 %20

# Opportunity costs 
sporeRep  <- mergedExpData_time_distinct %>%
  mutate(opportunity = 2*gene_length.filled*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(direct = 2*gene_length.filled*14)

# Sum 
sporeRepSum <- sporeRep %>%
  summarise(sumOpp = sum(opportunity, na.rm =T), sumDir = sum(direct, na.rm =T))
sporeRepTotal <- sporeRepSum$sumOpp+sporeRepSum$sumDir
# opportunity 50601460
# direct 20240584
# total 70842044

# percentage compared to total genome
sporeRepTotal/genome_tot  #17%

# Transcription costs
#####################
# Spore Transcription
#####################
##########################################################################

# 1 mRNA can yield to 100-1000 proteins (Cell biology by the numbers)
# I will count opportunity costs separately, so I can consider repolimerization costs

# Opportunity costs 
avg.protein.molecules <- 1774445
sporeTranscriptOpp <- mergedExpData_time_distinct %>%
  mutate(estimation = (((abundance.filled/1e2)*avg.protein.molecules)/1e6)*gene_length.filled) %>% #protein abundance/1000 X 1.8 X gene length
  # the reason I multiply with 1.8 is that the protein abundance is reported as parts per million. An average size bacteria has about 3 million protein molecules
  # this was reported experimentally as average 1774445 in Bacillus (Maass et.al. 2011)
  mutate(opportunity = estimation*31) 

# Sum 
sporeTranscriptOppSum <- sporeTranscriptOpp %>%
  group_by(time) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T))

# Sum stages
sporeTranscriptOppSum2 <- sporeTranscriptOpp %>%
  group_by(sourceID) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T))
  
# Direct costs 
sporeTranscriptDir <- mergedExpData_time %>%
  mutate(estimation = (((abundance.filled/1e2)*avg.protein.molecules)/1e6)*gene_length.filled) %>% 
  mutate(direct = estimation*(10+(2*12*1)))#hours
  # average sporulation time is 8hours, median mRNA degradation rate of Bacillus is 12 per hour 
  # (DOI: 10.1007/s00438-003-0883-6)
  # = 12 re-polymerization events per hour
  # assuming nucleotides are well recycled and it only affects polymerization costs

# Sum 
sporeTranscriptDirSum <-sporeTranscriptDir %>% 
  group_by(time_h) %>%
  summarise(sumDir = sum(direct, na.rm =T))

# Stages 
sporeTranscriptDirSum2 <-sporeTranscriptDir %>% 
  group_by(sourceID) %>%
  summarise(sumDir = sum(direct, na.rm =T))

# Cumulative costs
transcriptCosts <- cbind.data.frame(sporeTranscriptDirSum$time_h, opportunity = sporeTranscriptOppSum$sumOpp, 
                                    direct = sporeTranscriptDirSum$sumDir, 
                                    total = sporeTranscriptOppSum$sumOpp + sporeTranscriptDirSum$sumDir)

transcriptSum <- colSums(transcriptCosts[,-1])

#Cost of genes 
sporeTranscriptDirDist <- sporeTranscriptDir %>%
              group_by(locus_tag) %>%
              summarise(sumDist = sum(direct, na.rm =T))

# Translation costs
#####################
# Spore Translation
#####################

# Fill missing protein sequence cost estimations
median.aa_opportunitySum <- median(mergedExpData_time_distinct$aa_opportunitySum, na.rm = T) #5723
median.aa_directSum <- median(mergedExpData_time_distinct$aa_directSum, na.rm = T) #1351

# Opportunity and direct costs 
sporeTranslationOppDir <- mergedExpData_time_distinct %>%
  mutate(estimation = (abundance.filled*avg.protein.molecules)/1e6) %>% 
  mutate(aa_opportunitySum.filled = replace_na(aa_opportunitySum, median.aa_opportunitySum)) %>%
  mutate(aa_directSum.filled = replace_na(aa_directSum, median.aa_directSum)) %>%
  mutate(direct = estimation*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
  mutate(total = direct + opportunity)

# Sum
sporeTranslationOppDirSum <- sporeTranslationOppDir %>%
  group_by(time_h) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T), 
            total = sum(total, na.rm = T))

translationSum <- colSums(sporeTranslationOppDirSum[,-1])

# Opportunity 
opp <- translationSum[1]+transcriptSum[1]+septumOpp+genome_opp
# Direct
dir <- translationSum[2]+transcriptSum[2]+septumDir+genome_dir

(opp/(opp+dir))*100 #66.4
(dir/(opp+dir))*100 #33.6

# Total costs, plots, pie, bars
# Following will yield Figure 1
#####################
# Total, model, plot 
#####################

# Total costs 
cost_rep_all    <- genome_opp+genome_dir
cost_rep_part   <- sporeRepSum$sumOpp+sporeRepSum$sumDir
cost_rep_rest   <- cost_rep_all-cost_rep_part
cost_transcript <- transcriptSum[3]
cost_translation<- translationSum[3]
cost_membrane   <- septumOpp+septumDir
all_pie_costs   <- cost_rep_all+cost_transcript+cost_translation+cost_membrane

part_costs <- cost_rep_part+cost_transcript+cost_translation
# % proportions 
rep_partial     <- cost_rep_part/all_pie_costs*100 
rep_rest        <- cost_rep_rest/all_pie_costs*100 
rep_all         <- ((cost_rep_all)/all_pie_costs)*100
transcript      <- cost_transcript/all_pie_costs*100
translation     <- cost_translation/all_pie_costs*100
membrane        <- cost_membrane/all_pie_costs*100
  
proportion <- c(rep_all, transcript, translation, membrane)
pieCost    <- c("replication", "transcription", "translation", "membrane")
pieData    <- cbind.data.frame(pieCost, proportion)

pieData$labels <- paste(pieData$pieCost, round(pieData$proportion, 1), "%")

pieSpore <- ggplot(pieData, aes(x = "", y = proportion, fill = pieCost))+
  geom_bar(width = 1, stat = "identity")+ 
  coord_polar("y", start = 0)+
  mytheme+ 
  scale_fill_manual(values = c("#CC79A7", "#009E73","#D55E00","#0072B2"))+
  geom_text_repel(aes(label = labels), size = 4.5, show.legend = FALSE)+
  theme_void()+
  theme(legend.position = "none")

#ggsave("figures/figure1_pieSpore.pdf", pieSpore, height = 5, width = 6)

### Total costs of sporulation ###
### Figure  ###
time        <- rep(c(1:8), times = 2)
opportunity <- transcriptCosts$opportunity + sporeTranslationOppDirSum$opportunity
direct      <- transcriptCosts$direct + sporeTranslationOppDirSum$direct
costs       <- c(opportunity, direct)
type        <- rep(c("opportunity", "direct"), each = 8) 

sporulationCosts <- cbind.data.frame(time, costs, type)
sum(sporulationCosts$costs) #1554146274 (transcription and translation)

cost_rep_part+cost_rep_rest+sum(sporulationCosts$costs) #1967197654

sporulationTotal <- cbind.data.frame(time = c(1:8), 
                                     costs = direct + opportunity)

# Figure 1 Model 
fit <- nls(costs ~ SSasymp(time, yf, y0, log_alpha), data = sporulationTotal)
coef(fit) 

label    <- coef(fit)[3]

tt       <- seq(1,8, by = 0.1)
pred     <- predict(fit, list(time = tt))
preddata <- cbind.data.frame(pred, tt)

RSS.p <- sum(residuals(fit)^2)
y     <- as.numeric(as.character(sporulationTotal$costs))
TSS   <- sum((log10(y) - mean(y))^2)
R2    <- 1 - (RSS.p/TSS)

#scientific_10 <- function(x) {
#  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
#}

sporulationCosts_lay1 <- sporulationCosts %>%
  group_by(time) %>%
  summarize(sum = sum(costs))

sporulationCosts$type <- factor(sporulationCosts$type, c("opportunity","direct"))
# Plot 
# Note that some adjustments are done in Adobe Illustrator

my_lab <- c(expression(P['D']),
            expression(P['O']), 
            expression(P['T']))


f1 <- ggplot(NULL, aes(x = x, y = y))+
  geom_vline(xintercept = 2, linetype = "dashed")+
  geom_bar(data = sporulationCosts_lay1, 
           aes(x = time, y = sum), stat = "identity", color = "grey90", fill = "grey75", alpha = 0.5)+
  geom_bar(data = sporulationCosts, 
           aes(x = time, fill = type, y = costs), stat = "identity", color = "grey25", position = position_dodge(width=1))+
  ylab("ATP molecules")+
  xlab("Time (h)")+
  geom_line(data = preddata, aes(x = tt, y = pred), linewidth = 1)+
  mytheme+
  scale_y_continuous(breaks = c(2*10^8, 4*10^8, 6*10^8, 8*10^8, 10*10^8), 
                     labels = c(2,4,6,8,10), sec.axis=dup_axis())+
  scale_x_continuous(sec.axis=dup_axis())+
  annotate("text",x=-0.7,y=1.2*10^9,label=paste("(x10^8)"), parse =T, size = 18/.pt)+
  coord_cartesian(xlim = c(0.5, 8.5), clip="off")+
  
  annotate(geom = "text", x = 8, y = 2e8, label = paste("-\U03BB==", round(label, 3)), hjust = "right", size = 6, fontface = 'italic', parse = T)+
  annotate(geom = "text", x = 8, y = 1.2e8, label = paste("R^2==", round(R2, 3)), hjust = "right", size = 6, fontface = 'italic', 
           parse=TRUE)+
  theme(legend.position = c(0.32, 0.83), legend.title = element_blank())+
  scale_fill_manual(values = c("#D55E00","#0072B2"), labels=c(my_lab[1], 
                               my_lab[2],
                               my_lab[3]))

# ggsave("figures/figure1_sporeCostsTime.pdf", f1, height = 5, width = 6)
# END Figure 1 

# Costs until full commitment 
line_integral <- function(x, y) {
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx *my)
} 

x <- preddata$tt[1:11] #2 hours
y <- preddata$pred[1:11]
plot(x,y,"l")
commitment <- line_integral(x,y)
(commitment/sum(sporulationTotal))*100 #34.9%. %only transcription and translation 

   
                          ####################
########################## SPORE REVIVAL COSTS ##############################
                          ####################

####################################################
# Newly synthesized proteins from Swarge et al. 2020
####################################################

germination6_interval <- cbind.data.frame(
  protID = germination6$ProtID, H0.25 = germination6$T15-germination6$T0,
  H0.5 = germination6$T30-germination6$T15, H1 = germination6$T60-germination6$T30, 
  H1.5 = germination6$T90-germination6$T60, H2.5 = germination6$T150-germination6$T90, 
  H3.5 = germination6$T210-germination6$T150, H5.5 = germination6$T330-germination6$T210)

germLong6_interval <- germination6_interval %>%
  pivot_longer(cols = c(H0.25:H5.5),
               names_to =  "time_interval", values_to = "score") %>%
  filter(!score <= 0)

protSeqTidyAbun$gene <- protSeqTidyAbun$gene.y

germLong6_merged_interval <- germLong6_interval %>%
  left_join(protSeqTidyAbun, by = "protID") %>%
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(score.filled = replace_na(score, 18)) %>%
  mutate(protein_length.filled = replace_na(protein_length, 254)) %>%
  mutate(gene_length.filled = replace_na(gene_length, 765)) %>%
  mutate(time_h = case_when(
    time_interval == "H0.25" ~ 0.25,
    time_interval == "H0.5" ~ 0.25,
    time_interval == "H1" ~ 0.5,
    time_interval == "H1.5" ~ 0.5,
    time_interval == "H2.5" ~ 1,
    time_interval == "H3.5" ~ 1,
    time_interval == "H5.5" ~ 2))

### Replication costs ###

germRep6  <- germLong6_merged_interval %>%
  mutate(opportunity = 2*gene_length.filled*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(direct = 2*gene_length.filled*14)%>%
  distinct(protID, .keep_all = T) # keping distinct proteins not to account them multiple times 

# Sum 
germRepSum6 <- germRep6 %>%
  group_by(time_interval) %>%
  summarise(sumOpp = sum(opportunity, na.rm =T), sumDir = sum(direct, na.rm =T))

germRepTotal <- sum(germRepSum6$sumOpp+germRepSum6$sumDir)
  
germinationRepCost <- germRepSum6$sumOpp[1] + germRepSum6$sumDir[1]
outgrowthRepCost   <- sum(germRepSum6$sumOpp[2:7]) + sum(germRepSum6$sumDir[2:7])  
  
germRepTotal/genome_tot

### Transcription costs ###
# Opportunity costs 
germination6opp_interval   <- germLong6_merged_interval %>%
  mutate(estimation  = (score.filled/1e2)*(avg.protein.molecules)/1e6) %>% 
  mutate(opportunity = estimation*31) %>%
  group_by(time_interval) %>% 
  summarise(value  = sum(opportunity, na.rm = T))%>%
  mutate(source = rep("transcription"))%>%
  mutate(name = rep("opportunity")) 

# Direct costs 
germination6dir_interval    <- germLong6_merged_interval  %>%
  mutate(estimation  = (score.filled/1e2)*(avg.protein.molecules/1e6)*gene_length.filled) %>% 
  mutate(direct      = estimation*(10+(2*12*time_h)))%>%
  group_by(time_interval) %>% 
  summarise(value  = sum(direct, na.rm = T)) %>%
  mutate(source = rep("transcription")) %>%
  mutate(name = rep("direct")) 

totalGermtranscript <- sum(germination6opp_interval$value)+sum(germination6dir_interval$value)

### Translation costs ###

# Opportunity and direct costs 
germination6translation_interval <- germLong6_merged_interval  %>%
  mutate(estimation = ((score.filled)*(avg.protein.molecules))/1e6) %>% 
  mutate(direct = estimation*aa_directSum.filled) %>% # 
  mutate(opportunity = estimation*aa_opportunitySum.filled) %>% 
  mutate(total = direct + opportunity)%>% 
  group_by(time_interval) %>%
  summarise(opportunity = sum(opportunity, na.rm = T), 
            direct = sum(direct, na.rm = T)) %>%
  pivot_longer(cols = 2:3) %>%
  mutate(source = rep("translation"))

totalGermtranslation <- sum(germination6translation_interval$value)

germination6_all_interval <- rbind.data.frame(germination6dir_interval, germination6opp_interval, germination6translation_interval)

germination6sum_interval <- germination6_all_interval %>%
  group_by(time_interval, name) %>%
  summarise(costs = sum(value)) %>%
  mutate(hour = gsub("[A-Z]", "", time_interval)) 

germination6sum_interva_all <- germination6_all_interval %>%
  group_by(source) %>%
  summarise(costs = sum(value))

germinationTT <- sum(germination6sum_interval$costs[1:2])
outgrowthTT   <- sum(germination6sum_interval$costs[3:14])
outTT <- sum(germination6sum_interva_all$costs)+membraneGerm #1.4x10^10
germinationTT/(germinationTT+outgrowthTT)
germinationTT/outTT

germinationTT/(4321265075+germinationTT)

# Following will yield Figure 2
####################
# Model, plot, pie 
####################

germination6sum_interval$hours <- as.numeric(germination6sum_interval$hour)

#Add germination to the first data point 
germination6sum_interval[3,3] <- germination6sum_interval[3,3] + germination6sum_interval[1,3]
germination6sum_interval[4,3] <- germination6sum_interval[4,3] + germination6sum_interval[2,3]

germination6sum_int_sum <- germination6sum_interval %>%
  group_by(hours) %>%
  summarize(costs = sum(costs)) %>%
  filter(!hours == 0.25) %>%
  filter(!hours == 5.5)

fit3 <- nls(costs ~ SSasymp(hours, yf, y0, log_alpha), data = germination6sum_int_sum)
coef(fit3) 

label3    <- coef(fit3)[3]

tt3       <- seq(0.5, 4.9, by = 0.01)
pred3     <- predict(fit3, list(hours = tt3))
preddata3 <- cbind.data.frame(pred3, tt3)

RSS.p3 <- sum(residuals(fit3)^2)
y3     <- as.numeric(as.character(germination6sum_int_sum$costs))
TSS3   <- sum((log10(y3) - mean(y3))^2)
R23    <- 1 - (RSS.p3/TSS3)

# Plot 

germination6sum_interval_plotSum <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  filter(!hours == 5.5) %>%
  group_by(hours) %>%
  summarize(sum = sum(costs))

germination6sum_interval$name <- factor(germination6sum_interval$name, c("opportunity", "direct"))
  
germination6sum_interval_plotSum_n <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  filter(!hours == 5.5)

germination6sum_interval_plotSum_out <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  group_by(hours) %>%
  summarize(sum = sum(costs))

germination6sum_interval_plotSum_out_all <- germination6sum_interval %>%
  filter(!hours == 0.25) %>%
  group_by(hours, name) %>%
  summarize(sum = sum(costs))

min15_sep <- germination6sum_interval %>%
  filter(hours == 0.25)

# Figure 1 - germination - new 
germs <- ggplot(NULL, aes(x = x, y = y))+
  geom_vline(xintercept = 4.82, linetype = "dashed")+
  geom_bar(data = germination6sum_interval_plotSum_out, 
           aes(x = seq(0.5, by = mean(diff(hours)), length = length(hours)), 
                                   y = sum), stat = "identity", color = "grey90", fill = "grey75", alpha = .5)+
  
  annotate("rect", xmin = 0.11, xmax = 0.85, ymin = 0, ymax = 84970721+323997094,
           alpha = .5,fill = "grey25")+
  geom_segment(aes(x = 0.11, xend = 0.85, y = 84970721+323997094, yend = 84970721+323997094), color = "grey25")+
  geom_bar(data = germination6sum_interval_plotSum_out_all, 
           aes(x = seq(0.5, by = mean(diff(hours)), length = length(hours)), fill = name, 
               y = sum), position = position_dodge(0.9), stat = "identity", color = "grey25")+

  ylab("ATP molecules")+
  xlab("Time (h)")+
  geom_line(data = preddata3, aes(x = tt3, y = pred3), linewidth = 1)+
  annotate(geom = "text", x = 4.8, y = 5.1e8, label = paste("-\U03BB==", round(label3, 3)), hjust = "right", size = 5, fontface = 'italic', parse = T)+
  annotate(geom = "text", x = 4.8, y = 3.8e8, label = paste("R^2==", round(R23, 3)), hjust = "right", size = 5, fontface = 'italic', 
           parse=TRUE)+
  mytheme+
  scale_y_continuous(breaks = c(4*10^8, 8*10^8, 12*10^8, 16*10^8, 20*10^8), 
                     labels = c(4,8,12,16,20), sec.axis=dup_axis())+
  scale_x_continuous(sec.axis=dup_axis())+
  annotate("text",x=-0.55,y=2.3*10^9,label=paste("(x10^8)"), parse =T, size = 16/.pt)+
  coord_cartesian(xlim = c(0.25, 6), clip="off")+
  scale_fill_manual(values = c("#D55E00","#0072B2"), labels=c(my_lab[1], 
                          my_lab[2],
                          my_lab[3]))+
  annotate("rect", xmin = 0.33, xmax = 0.67 , ymin = 0, ymax = 323997094,
           alpha = .6,fill = "grey25")+
  geom_segment(aes(x = 0.33, xend = 0.67, y = 323997094, yend = 323997094), color = "grey25")+
  annotate("rect", xmin = 0.71, xmax = 1.05, ymin = 0, ymax = 84970721,
                                              alpha = .6,fill = "grey25")+
  geom_segment(aes(x = 0.71, xend = 1.05, y = 84970721, yend = 84970721), color = "grey25")+
  theme(legend.position = c(0.22, 0.88), legend.title = element_blank())

#ggsave("figures/figure2_germCostsTime.pdf", germs, height = 5, width = 6)
# End

# Pie
replicationGerm   = germinationRepCost+outgrowthRepCost
transcriptionGerm = totalGermtranscript
translationGerm   = totalGermtranslation
membrGerm         = membraneGerm

allGerm <- transcriptionGerm+translationGerm+membrGerm

transcriptionGerm = totalGermtranscript/outTT*100 
translationGerm   = totalGermtranslation/outTT*100 
membrGerm         = membraneGerm/outTT*100 #%67.9

proportion2 <- c(round(transcriptionGerm, 2), round(translationGerm, 2), round(membrGerm, 2))
pieCost2    <- c("transcription", "translation", "membrane")
pieData2    <- cbind.data.frame(pieCost2, proportion2)

pieData2$labels <- paste(pieData2$pieCost, round(pieData2$proportion, 1), "%")

pieGerm <- ggplot(pieData2, aes(x = "", y = proportion2, fill = pieCost2))+
  geom_bar(width = 1, stat = "identity")+ 
  coord_polar("y", start=0)+
  mytheme+ 
  scale_fill_manual(values = c("#CC79A7","#D55E00","#0072B2"))+
geom_text_repel(aes(label = labels), size = 4.5, show.legend = FALSE)+
theme_void()+
theme(legend.position = "none")

#ggsave("figures/figure2_pieGerm.pdf", pieGerm, height = 5, width = 5)

# Following will yield Figure 3
####################################
# HEAD TO HEAD COMPARISON OF TRAITS 
####################################

# Merge trait data
mergedTraitData <- otherTraits %>%
  filter(!category == "sporulation") %>%
  filter(!category == "germination") %>%
  left_join(annotationData, by = "gene") %>%
  left_join(protSeqTidyAbun[,c("locus_tag","abundance", "aa_opportunitySum", "aa_directSum")], by = "locus_tag") %>%
  mutate(abundance.filled = replace_na(abundance, 18)) %>%
  group_by(category) %>%
  mutate(no_genes = length(category))

# All costs
totalCosts_traits <- mergedTraitData %>%
  mutate(translationAll = abundance.filled*(avg.protein.molecules/1e6)) %>%
  mutate(aa_opportunitySum.filled = median(aa_opportunitySum, na.rm = T)) %>%
  mutate(aa_directSum.filled = median(aa_directSum, na.rm = T)) %>%
  mutate(translationDirect = translationAll*aa_directSum.filled) %>% # ignoring protein degradation
  mutate(translationOpportunity = translationAll*aa_opportunitySum.filled) %>%
  mutate(translationTotal = translationDirect + translationOpportunity) %>%
  mutate(transcriptionAll = (abundance.filled/1e2)*(avg.protein.molecules/1e6)*as.numeric(gene_length)) %>%
  mutate(transcriptionDirect = transcriptionAll*(10+(2*12*1))) %>% #assuming that mRNAs transcribed at least 1 hour
  mutate(transcriptionOpportunity = transcriptionAll*31) %>%
  mutate(transcriptionTotal = transcriptionDirect + transcriptionOpportunity) %>%
  mutate(RepOpportunity = 2*as.numeric(gene_length)*35) %>% #2 = doublestring  35 = nucleotide costs 
  mutate(RepDirect = 2*as.numeric(gene_length)*14) %>%
  mutate(RepTotal = RepOpportunity+RepDirect) %>%
  mutate(costs = transcriptionTotal + translationTotal+ RepTotal) # not accounting replicat

germinationTTrep = germinationTT+germinationRepCost
outgrowthTTrep   = outgrowthTT+outgrowthRepCost

# Cumulative costs
costs_sum_traits <- totalCosts_traits %>%
  group_by(category) %>%
  summarize(sumCosts = sum(costs, na.rm = T),
            number_genes = mean(no_genes)) %>%
  select(category, sumCosts) %>%
  add_row(category = "growth requirements", sumCosts = 2.6e10) %>%
  add_row(category = "basal metabolism", sumCosts = 3.5e8) %>%
  add_row(category = "membrane", sumCosts = membraneOpp+membraneDir) %>%
  add_row(category = "developmental programm", sumCosts = part_costs+
            germinationTTrep + outgrowthTTrep) %>%
  add_row(category = "genome replication", sumCosts = genome_tot) 

category = factor(c("sporulation", "germination", "outgrowth"))
sumCosts = c(part_costs,germinationTTrep,outgrowthTTrep)
costs_sum_dev <- cbind.data.frame(category, sumCosts)

# Add a second axis 
costs_sum_traits_rel <- costs_sum_traits %>%
  mutate(relative = (sumCosts/2.6e10)*100)

### Figure 3 ###
sumAxis = part_costs+germinationTTrep+outgrowthTTrep

part_costs/sumAxis*100
germinationTTrep/sumAxis*100
outgrowthTTrep/sumAxis*100

log10(sumAxis)
log10(outgrowthTTrep+germinationTTrep)

costs_sum_traits_rel$category2 <- c("biofilm structure", "chemotaxis", "competence", 
                                    "essential genes", "flagella", "heat-shock proteins", 
                                    "homeostasis", "swarming", "total cell budget", "basal metabolism h-1", 
                                    "membrane lipid synthesis", "developmental programm", "genome replication")
  
 f3 <- ggplot(costs_sum_traits_rel, aes(x = log10(sumCosts), y = reorder(category2, sumCosts)))+
  geom_col(fill="grey75", color="grey25")+
  geom_vline(xintercept = log10(germinationTTrep))+
  geom_vline(xintercept = log10(germinationTTrep+part_costs))+
  geom_vline(xintercept = log10(germinationTTrep+outgrowthTTrep+part_costs))+
  mytheme2+
  # Custom the Y scales:
  scale_x_continuous(
                  
    # Features of the first axis
    name = "ATP molecules",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(trans = ~log10((10^./2.6e10)*100), name = "Costs relative to the cell budget (%) "))+
  
  coord_cartesian(xlim = c(7.8,10.5))+
  
  theme(axis.ticks.y = element_blank())+
  theme(axis.title.y = element_blank())

# ggsave("figures/figure3.pdf", f3, height = 5, width = 7)

# germination, sporulation, outgrowth 
ggplot(costs_sum_dev, aes(y = log10(sumCosts), x = reorder(category, sumCosts)))+
  geom_col(fill="grey", color="black")+
  scale_y_continuous(name = "Costs (ATP molecules)")+
  coord_cartesian(ylim = c(7.8,10.5))+
  mytheme

