# COGs and evolution plots
# 26 August 2025 - last update
# Author: C. Karakoc

######################
# Packages & Plotting 
######################

library(tidyverse)       
library(stringr) 
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

######################
# Data
######################

# COGs - Galperin et al., 2022
meta_data   <- read.table("./evolution/empirical/data/calc_meta.csv", sep = ",", dec = "." , header = T, stringsAsFactors = F)

meta_data$No._proteins_numeric <- as.numeric(gsub(",", "", meta_data$No._proteins))

meta_data_clean <- meta_data %>% 
  filter(!Spore_former == "Unk") %>% 
  mutate(allNorm = All_spore_COGs_.out_of_237./ No._proteins_numeric,
         narrowNorm = Narrowly_conserved_.out_of_74. / No._proteins_numeric,
         coreNorm = Core_COGs_.out_of_112./ No._proteins_numeric, 
         allNorm_genome = All_spore_COGs_.out_of_237./ Genome_size_Mb,
         narrowNorm_genome = Narrowly_conserved_.out_of_74. / Genome_size_Mb,
         coreNorm_genome = Core_COGs_.out_of_112./ Genome_size_Mb)

msPlot <- ggplot(meta_data_clean, aes(x = All_spore_COGs_.out_of_237., color = Spore_former))+
  geom_histogram(aes(y = ..density.., fill = Spore_former),color = "grey50", bins = 30, alpha = 0.4, linewidth = 0.5)+
  geom_density(linewidth = 1.2, aes(group = Spore_former))+
  xlab("Present spore COGs")+
  ylab("Frequency")+
  mytheme+
  scale_fill_manual(values = c("#0072B2", "#D55E00"), labels = c("Non-spore-former", "Spore-former")) +
  scale_color_manual(values = c("#0072B2", "#D55E00"), labels = c("Non-spore-former", "Spore-former")) +
  theme(legend.position = c(0.75,0.85))+
  theme(legend.text = element_text(size = 14))+
  scale_y_continuous(sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())
                    
#ggsave("./evolution/empirical/figures/COGS.pdf", plot = msPlot, width = 5.6, height = 4.4, dpi = 300)

############################################
# Plot's for Model results for the main text
############################################

evoratio <- read.table("./evolution/model/data/evo_ratio_conditional.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
evoratio$values <- strsplit(evoratio$V2, ",")

max_values_evo <- max(sapply(evoratio$values, length))
value_cols_evo <- paste0("value", 1:max_values_evo)
values_df_evo  <- as.data.frame(do.call(rbind, lapply(evoratio$values, `length<-`, max_values_evo)))
colnames(values_df_evo) <- value_cols_evo

new_data_evo <- cbind(evoratio[,1, drop=FALSE], values_df_evo)

new_data_matrix_evo <- as.matrix(new_data_evo[, -1])
rownames(new_data_matrix_evo) <- new_data_evo[, 1]
new_data_transposed_evo <- t(new_data_matrix_evo)

new_data_final_evo <- as.data.frame(new_data_transposed_evo, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final_evo)) {
  new_data_final_evo[, i] <- as.numeric(new_data_final_evo[, i])
}

new_data_final_evo$average_Ne <- (new_data_final_evo$ratio_rates_lower_Ne+new_data_final_evo$ratio_rates_higher_Ne)/2

# Custom labels using expression()
custom_labels <- c(
  expression(italic("Neutrality")),
  expression(italic(N[e]) ~ "= 3.2 × 10"^8),
  expression(italic(N[e]) ~ "= 6.1 × 10"^8)
)

evorates_av <- new_data_final_evo %>%
  select(deletion_size, null_ratio_rates, average_Ne) %>%
  pivot_longer(-deletion_size, names_to = "rates", values_to = "values") 

custom_labels2 <- c(
  expression("Average empirical " * N[e]),
  "Neutrality"
)

evorates_av_plot <- ggplot(evorates_av, aes(x = log10(deletion_size), y = log10(values), fill = rates, color = rates))+
  geom_hline(yintercept = 0, linetype = "dotted", size = 1, color = "grey25")+
  geom_point(size = 3, shape = 21, color = "grey25")+
  geom_line(size = 1.2, aes(linetype = rates))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  xlab(expression("Deletion size (bp),"~Delta))+
  ylab(expression(atop("Ratio of deletion size " * Delta, "relative to substitutions, " * frac(d[del](Delta), d[sub]))))+
  mytheme+
  theme(legend.position = c(0.3, 0.85), legend.background = element_rect(fill = "white"))+
  theme(legend.text = element_text(size = 14))+
  scale_fill_manual(values = c("#0072B2", "grey25"), labels = custom_labels2)+
  scale_color_manual(values = c("#0072B2", "grey25"), labels = custom_labels2)+
  scale_x_continuous(labels = function(x) parse(text = paste0("10^", format(x, digits = 1))),sec.axis = dup_axis())+
  scale_y_continuous(labels = function(x) parse(text = paste0("10^", format(x, digits = 1))),sec.axis = dup_axis())+
  guides(linetype = "none")

### Combine plots with empirical data ###
# These will yield Figure 5 

figure5 <- ggarrange(msPlot, evorates_av_plot, labels = c("A", "B"),
                     ncol = 2, nrow = 1)

#ggsave("./evolution/empirical/figures/figure5.pdf", plot = figure5, width = 13, height = 5, dpi = 300)
