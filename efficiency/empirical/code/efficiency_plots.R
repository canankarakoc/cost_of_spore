# Spore Costs
#26 August 2025 - last update
# Author: C. Karakoc

######################
# Packages & Plotting 
######################

library(tidyverse)       
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

setwd("~/Documents/GitHub/cost_of_spore")

######################
# Data
######################

# Sporulation efficiency (corresponds to the SI Table 1 with full dataset of the papers)
efficiencyEmp     <- read.table("./efficiency/empirical/data/efficiency.csv", sep = ",", dec = "," , header = T, stringsAsFactors = F)

#################################
# Empirical spore efficiency
#################################

efficiencyEmp$efficiency <- as.numeric(efficiencyEmp$efficiency)
med.efficiency <- median(as.numeric(efficiencyEmp$efficiency))

efficiencyPlot <- ggplot(efficiencyEmp, aes(efficiency))+
  geom_density(linewidth = 1.6, color = "#0072B2")+
  xlab("Sporulation efficiency (%)")+
  ylab("Frequency")+
  geom_vline(xintercept = med.efficiency, color = "black", linetype = "dashed" )+
  coord_cartesian(ylim = c(0.0025, 0.0125))+
  annotate(geom = "label", x = 33, y = 0.005, color ="black", fill = "white", label = "Median = 31%")+
  scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,75,100), labels = c(0,25,50,75,100))+
  mytheme
  
library(truncnorm)
fitE   <- density(efficiencyEmp$efficiency, from = 0, to = 100)
N     <- 1e6
x.new <- rtruncnorm(a = 0, b = 100, N, sample(efficiencyEmp$efficiency, size = N, replace = TRUE))
plot(density(x.new, bw = fitE$bw))

densityPlot <- as.data.frame(x.new)

means = densityPlot %>%
  summarise(M = median(x.new), SD = sd(x.new), N = n())

# Sample size of 50. T Distribution intervals
means$error = qt(0.975, df = 50-1)*means$SD/sqrt(means$N)
means$upper = means$M+means$error
means$lower = means$M-means$error

efficiencyPlot <- ggplot(densityPlot, aes(x = x.new))+
  geom_density(alpha = 0.1, linewidth = 0.8, bw = 11.35555, color = "#0072B2")+
  geom_vline(xintercept = med.efficiency, linetype = 'dashed', color = "black")+
  annotate(geom = "label", x = 33, y = 0.005, color = "black", fill = "white", label = "Median = 31%")+
  labs(x="Sporulation efficiency", y = "Density")+
  mytheme + 
  scale_x_continuous(limits = c(-30,130), sec.axis = dup_axis())+
  scale_y_continuous(sec.axis = dup_axis())

#ggsave("./efficiency/empirical/figures/efficiencyEmp.pdf", efficiencyPlot, width = 5.6, height = 4.4)

############################################
# Plot's for Model results for the main text
############################################

# Consumer-resource model 
consumerResource <- read.table("./efficiency/model/data/consumer_resource_spore_model.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
consumerResource$values <- strsplit(consumerResource$V2, ",")

max_values <- max(sapply(consumerResource$values, length))
value_cols <- paste0("value", 1:max_values)
values_df  <- as.data.frame(do.call(rbind, lapply(consumerResource$values, `length<-`, max_values)))
colnames(values_df) <- value_cols

new_data <- cbind(consumerResource[,1, drop=FALSE], values_df)

new_data_matrix <- as.matrix(new_data[, -1])
rownames(new_data_matrix) <- new_data[, 1]
new_data_transposed <- t(new_data_matrix)

new_data_final <- as.data.frame(new_data_transposed, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final)) {
  new_data_final[, i] <- as.numeric(new_data_final[, i])
}

legend_names <- c("resource concentration", "cell density", "spore density")  # Replace with your desired legend names

plotModel <- new_data_final%>%
  pivot_longer(-hours, names_to = "population", values_to = "abundance")

plotModel$population <- factor(plotModel$population, levels = c("resource_concentration", "cell_concentration", "spore_concentration"))

annotation_positions <- data.frame(
  population = c("resource", "vegetative", "spore"),
  hours = c(1, 30, 34),  # Example positions for each level
  abundance = c(10^3.5, 10^9.5, 10^8.5)  # Example positions for each level
)

plot1 <- ggplot(plotModel, aes(x = hours, y = abundance, color = population))+
  geom_line(size = 1.6)+
  mytheme+
  labs(x = "Time (h)", y = "Concentration")+
  theme(legend.position = c(0.25, 0.88))+
  scale_color_manual(labels = c("resource concentration", "cell density", "spore density"), 
                     values = c("#009E73", "#0072B2", "#D55E00"))+
  scale_y_continuous(limits = c(10^0, 10^10), trans = "log10", breaks = c(10^1,10^3,10^5,10^7,10^9), 
                     labels = function(x) parse(text = paste0("10^", format(log10(x), digits = 2))), sec.axis = dup_axis())+
  scale_x_continuous(sec.axis = dup_axis())+
  theme(legend.position = "none")+
  geom_text(data = annotation_positions,
            aes(label = population), color = c("#009E73", "#0072B2", "#D55E00"),
            x = annotation_positions$hours, y = log10(annotation_positions$abundance),
            hjust = -0.1, vjust = 0.5, size = 6)   

#############################################################################

# Spore efficiency 
efficiency <- read.table("./efficiency/model/data/efficiency_data.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
efficiency$values <- strsplit(efficiency$V2, ",")

max_values_eff <- max(sapply(efficiency$values, length))
value_cols_eff <- paste0("value", 1:max_values_eff)
values_df_eff  <- as.data.frame(do.call(rbind, lapply(efficiency$values, `length<-`, max_values_eff)))
colnames(values_df_eff) <- value_cols_eff

new_data_eff <- cbind(efficiency[,1, drop=FALSE], values_df_eff)

new_data_matrix_eff <- as.matrix(new_data_eff[, -1])
rownames(new_data_matrix_eff) <- new_data_eff[, 1]
new_data_transposed_eff <- t(new_data_matrix_eff)

new_data_final_eff <- as.data.frame(new_data_transposed_eff, stringsAsFactors = FALSE)
for (i in 1:ncol(new_data_final_eff)) {
  new_data_final_eff[, i] <- as.numeric(new_data_final_eff[, i])
}

# spore/cell empirical estimate
ratio_estimate <- 2385796113/93854440000 # from bioaccounting C_spore/C_cell_budget

efficiency_estimate <- ggplot(new_data_final_eff,
            aes(x = spore_vs_cell_atp_ratio, y = sporulation_efficiency)) +
  geom_line(size = 1.3, colour = "#0072B2") +
  geom_vline(xintercept = ratio_estimate, linetype = "dashed") +
  labs(x = "Ratio of energetic costs, spore/cell",
       y = "Sporulation efficiency") +
  scale_x_log10(breaks = 10^seq(-4, 4, 1), labels = label_log()) +
  scale_y_log10(limits = c(1e-3, 1e-1),
                breaks = 10^seq(-3, 0, 1),
                labels = label_log()) +
  mytheme+
  annotate("text", x = 10^-1.9, y = 10^-2.5, label = "Emprical estimate", size = 6, angle = 90)
##########################################################################################################

figure4 <- ggarrange(efficiencyPlot, plot1, efficiency_estimate,
                     labels = c("A", "B", "C"),
                     ncol = 3, nrow = 1)

#ggsave("./efficiency/empirical/figures/figure4.pdf", plot = figure4, width = 15, height = 4, dpi = 300)

