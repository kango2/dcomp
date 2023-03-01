##Loading the data

##Loading the data
#Loading H4K20me1 (female) data
library(readr)
F_H4K20me1 <- read_delim("F_H4K20me1_peaks.broadPeak", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

colnames(F_H4K20me1) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value")

#Loading H4K20me1 (male) data
M_H4K20me1 <- read_delim("M_H4K20me1_peaks.broadPeak", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

colnames(M_H4K20me1) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value")

##Combining the data so it can be plotted together
#Making one dataframe with both samples
F_H4K20me1$Sample <- 'Female'
M_H4K20me1$Sample <- 'Male'
total_H4K20me1 <- rbind(F_H4K20me1, M_H4K20me1)

#Making new fold change column to manipulate
total_H4K20me1$Adjusted_fold_change <- total_H4K20me1$Fold_change

#Making female fold changes negative
total_H4K20me1[total_H4K20me1$Sample == 'Female', "Adjusted_fold_change"] <- -total_H4K20me1[total_H4K20me1$Sample == 'Female', "Adjusted_fold_change"]

#Making useful dataframe
Usable_H4K20me1 <- total_H4K20me1[c(-4, -5, -6, -8)]
Usable_H4K20me1 <- subset(Usable_H4K20me1, Usable_H4K20me1$`log10_Q-value` > 1.3)
#Filters for q value below 0.05
Usable_H4K20me1$Modification <- "H4K20me1"

#

#####GERMAN BELOW

##Set working directory
setwd("/Users/ashleymilton/OneDrive - UNSW/2021 PhD/Platypus ChIP/macs2_German/")


##Loading the data
#Loading H3K4me3 (female) data
library(readr)
F_H3K4me3 <- read_delim("K4me3_F_354_narrow_peaks.narrowPeak", 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)

colnames(F_H3K4me3) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value", "Peak")

#Loading H3K4me3 (male) data
M_H3K4me3 <- read_delim("K4me3_M_324_narrow_peaks.narrowPeak", 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)

colnames(M_H3K4me3) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value", "Peak")

##Combining the data so it can be plotted together
#Making one dataframe with both samples
F_H3K4me3$Sample <- 'Female'
M_H3K4me3$Sample <- 'Male'
total_H3K4me3 <- rbind(F_H3K4me3, M_H3K4me3)

#Making new fold change column to manipulate
total_H3K4me3$Adjusted_fold_change <- total_H3K4me3$Fold_change

#Making female fold changes negative
total_H3K4me3[total_H3K4me3$Sample == 'Female', "Adjusted_fold_change"] <- -total_H3K4me3[total_H3K4me3$Sample == 'Female', "Adjusted_fold_change"]

#Making useful dataframe
Usable_H3K4me3 <- total_H3K4me3[c(-4, -5, -6, -8, -10)]
Usable_H3K4me3 <- subset(Usable_H3K4me3, Usable_H3K4me3$`log10_Q-value` > 1.3)
#Filters for q value below 0.05
Usable_H3K4me3$Modification <- "H3K4me3"

#

##Loading the data
#Loading H3K27ac (female) data
library(readr)
F_H3K27ac <- read_delim("K27ac_F_330_narrow_peaks.narrowPeak", 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)

colnames(F_H3K27ac) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value", "Peak")

#Loading H3K27ac (male) data
M_H3K27ac <- read_delim("K27ac_M_300_narrow_peaks.narrowPeak", 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)

colnames(M_H3K27ac) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value", "Peak")

##Combining the data so it can be plotted together
#Making one dataframe with both samples
F_H3K27ac$Sample <- 'Female'
M_H3K27ac$Sample <- 'Male'
total_H3K27ac <- rbind(F_H3K27ac, M_H3K27ac)

#Making new fold change column to manipulate
total_H3K27ac$Adjusted_fold_change <- total_H3K27ac$Fold_change

#Making female fold changes negative
total_H3K27ac[total_H3K27ac$Sample == 'Female', "Adjusted_fold_change"] <- -total_H3K27ac[total_H3K27ac$Sample == 'Female', "Adjusted_fold_change"]

#Making useful dataframe
Usable_H3K27ac <- total_H3K27ac[c(-4, -5, -6, -8, -10)]
Usable_H3K27ac <- subset(Usable_H3K27ac, Usable_H3K27ac$`log10_Q-value` > 1.3)
#Filters for q value below 0.05
Usable_H3K27ac$Modification <- "H3K27ac"

#

##Loading the data
#Loading H3K4me (female) data
library(readr)
F_H3K4me <- read_delim("K4me_F_342_broad_peaks.broadPeak", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

colnames(F_H3K4me) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value")

#Loading H3K4me (male) data
M_H3K4me <- read_delim("K4me_M_312_broad_peaks.broadPeak", 
                         "\t", escape_double = FALSE, col_names = FALSE, 
                         trim_ws = TRUE)

colnames(M_H3K4me) <- c("Chromosome", "Start", "End", "ID", "Integer_score", "Strand", "Fold_change", "log10_P-value", "log10_Q-value")

##Combining the data so it can be plotted together
#Making one dataframe with both samples
F_H3K4me$Sample <- 'Female'
M_H3K4me$Sample <- 'Male'
total_H3K4me <- rbind(F_H3K4me, M_H3K4me)

#Making new fold change column to manipulate
total_H3K4me$Adjusted_fold_change <- total_H3K4me$Fold_change

#Making female fold changes negative
total_H3K4me[total_H3K4me$Sample == 'Female', "Adjusted_fold_change"] <- -total_H3K4me[total_H3K4me$Sample == 'Female', "Adjusted_fold_change"]

#Making useful dataframe
Usable_H3K4me <- total_H3K4me[c(-4, -5, -6, -8)]
Usable_H3K4me <- subset(Usable_H3K4me, Usable_H3K4me$`log10_Q-value` > 1.3)
#Filters for q value below 0.05
Usable_H3K4me$Modification <- "H3K4me"

#

#Making overall data frame for German

Platypus_ChIP_german <- rbind(Usable_H3K4me3, Usable_H3K27ac, Usable_H3K4me)

#Making overall old data frame

Platypus_ChIP <- rbind(Usable_H4K20me1)
Platypus_ChIP <- rbind(Platypus_ChIP, Platypus_ChIP_german)

#Log2 transformation
Platypus_ChIP$Log2_fold_change <- log2(Platypus_ChIP$Fold_change)

Platypus_ChIP$Log2_adjusted_fold_change <- Platypus_ChIP$Log2_fold_change
Platypus_ChIP[Platypus_ChIP$Sample == 'Female', "Log2_adjusted_fold_change"] <- -Platypus_ChIP[Platypus_ChIP$Sample == 'Female', "Log2_adjusted_fold_change"]





##Adding GFF
library(readr)
genes <- read_delim("~/Desktop/ChIP data local/Platypus_Genome/GCF_004115215.2_mOrnAna1.pri.v4_genomic.gff", 
                    "\t", escape_double = FALSE, col_names = FALSE, 
                    trim_ws = TRUE, skip = 8)


genes <- genes[which(genes[,3] == "gene"),]

library("reshape2")
genes[,9:14] <- colsplit(genes$X9, ";", c("1","2","3","4","5","6"))
genes[,15:16] <- colsplit(genes$"5", "=", c("7","Gene"))

genes <- genes[,c(1,4,5,7,16)]

colnames(genes) <- c("Chromosome", "Start", "End", "Strand", "Gene")

library(scales)
library(ggplot2)
library(gridExtra)
library(dplyr)

LISTOFCHROMOSOMES <- c("NC_041728.1", "NC_041729.1", "NC_041730.1", "NC_041731.1", "NC_041732.1", "NC_041733.1", "NC_041734.1", "NC_041735.1", "NC_041736.1", "NC_041737.1",
                       "NC_041738.1", "NC_041739.1", "NC_041740.1", "NC_041741.1", "NC_041742.1", "NC_041743.1", "NC_041744.1", "NC_041745.1", "NC_041746.1", "NC_041747.1",
                       "NC_041748.1", "NC_041749.1", "NC_041750.1", "NC_041751.1", "NC_041752.1", "NC_041753.1")

i <- "NC_041753.1"


chrom_sizes <- read_csv("/Users/ashleymilton/Desktop/ChIP data local/BigWig/ornAna1_chrom_sizes_NCBInames.csv", col_names = TRUE)
colnames(chrom_sizes) <- c("Chromosome", "Chr_Length")
Platypus_ChIP2 <- Platypus_ChIP %>% right_join(chrom_sizes, by=c("Chromosome"))
  

for(i in LISTOFCHROMOSOMES){
  genes_mod <- genes[which(genes$Chromosome == i),]
  Plot_data <- Platypus_ChIP2[which(Platypus_ChIP2$Chromosome == i),]
  chr_length <- unique(Plot_data$Chr_Length)
  plot <- ggplot(transform(Plot_data, Modification=factor(Modification, levels = c("H3K4me", "H3K27ac", "H3K4me3", "H4K20me1")))) +
    geom_rect(aes(xmin=Start,xmax=End,ymin=0,ymax=Log2_adjusted_fold_change,fill=Modification, colour=Modification)) +
    theme_classic() +
    #labs(caption = i)+
    scale_x_continuous(limits=c(0, chr_length), labels = label_number(scale = 1/1000000))+
    theme(plot.caption = element_text(hjust=0.5), legend.position = "none", strip.background = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), strip.text = element_blank()) +
    #theme(plot.caption = element_text(hjust=0.5), legend.position = "none", strip.background = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.line = element_blank(), axis.text = element_blank(),axis.ticks = element_blank(), strip.text.x = element_blank()) +
    scale_fill_manual(values = c("#58B38C50","#58B38C50","#58B38C50","#58B38C50")) +
    scale_color_manual(values = c("#58B38C50","#58B38C50","#58B38C50","#58B38C50")) +
    theme(panel.spacing.y = unit(-0.05, "cm")) +
    facet_wrap(~ Modification, ncol = 1, scales = "free_y")+
    geom_segment(data=genes_mod, aes(x=Start, xend=End, y=0, yend=0), colour="#6b6b6b", size=2)
  ggsave(plot=plot, filename = paste('yaxis_log2_facetfree_',i,'_q0.05_germaninput_noaxes.jpg'), path = "/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/ChIP_seq/Figures/Ash_ChIP_jpgs_030622_For_paper", width = 18, height = 6, dpi = 300)
}


###########################################

LISTOFHISTMODS <- c("H3K4me3", "H3K27ac", "H3K4me")


#Loop for each histone
for(i in LISTOFHISTMODS){
  x <- sum(Platypus_ChIP2$Modification == i & Platypus_ChIP2$Sample == 'Female')
  write.csv(x, file = paste(i,'_female_peaks.txt'))
  x <- sum(Platypus_ChIP2$Modification == i & Platypus_ChIP2$Sample == 'Male')
  write.csv(x, file = paste(i,'_male_peaks.txt'))
  x <- sum(Platypus_ChIP2$Modification == i & Platypus_ChIP2$Sample == 'Female' & (Platypus_ChIP2$Chromosome == 'NC_041749.1' | Platypus_ChIP2$Chromosome == 'NC_041750.1' | Platypus_ChIP2$Chromosome == 'NC_041751.1' | Platypus_ChIP2$Chromosome == 'NC_041752.1' | Platypus_ChIP2$Chromosome == 'NC_041753.1'))
  write.csv(x, file = paste(i,'_female_Xs_peaks.txt'))
  x <- sum(Platypus_ChIP2$Modification == i & Platypus_ChIP2$Sample == 'Male' & (Platypus_ChIP2$Chromosome == 'NC_041749.1' | Platypus_ChIP2$Chromosome == 'NC_041750.1' | Platypus_ChIP2$Chromosome == 'NC_041751.1' | Platypus_ChIP2$Chromosome == 'NC_041752.1' | Platypus_ChIP2$Chromosome == 'NC_041753.1'))
  write.csv(x, file = paste(i,'_male_Xs_peaks.txt'))
  x <- sum(Platypus_ChIP2$Modification == i & Platypus_ChIP2$Sample == 'Female' & Platypus_ChIP2$Chromosome != 'NC_041749.1' & Platypus_ChIP2$Chromosome != 'NC_041750.1' & Platypus_ChIP2$Chromosome != 'NC_041751.1' & Platypus_ChIP2$Chromosome != 'NC_041752.1' & Platypus_ChIP2$Chromosome != 'NC_041753.1')
  write.csv(x, file = paste(i,'_female_autos_peaks.txt'))
  x <- sum(Platypus_ChIP2$Modification == i & Platypus_ChIP2$Sample == 'Male' & Platypus_ChIP2$Chromosome != 'NC_041749.1' & Platypus_ChIP2$Chromosome != 'NC_041750.1' & Platypus_ChIP2$Chromosome != 'NC_041751.1' & Platypus_ChIP2$Chromosome != 'NC_041752.1' & Platypus_ChIP2$Chromosome != 'NC_041753.1')
  write.csv(x, file = paste(i,'_male_autos_peaks.txt'))
}


Platypus_ChIP_X1spec <- subset(Platypus_ChIP2, (Platypus_ChIP2$Chromosome == 'NC_041749.1' & Platypus_ChIP2$Start > 45100000))
Platypus_ChIP_X2spec <- subset(Platypus_ChIP2, (Platypus_ChIP2$Chromosome == 'NC_041750.1' & Platypus_ChIP2$Start > 5700000 & Platypus_ChIP2$End < 20000000))
Platypus_ChIP_X3spec <- subset(Platypus_ChIP2, (Platypus_ChIP2$Chromosome == 'NC_041751.1' & Platypus_ChIP2$Start > 18500000 & Platypus_ChIP2$End < 31900000))
Platypus_ChIP_X4spec <- subset(Platypus_ChIP2, (Platypus_ChIP2$Chromosome == 'NC_041752.1' & Platypus_ChIP2$Start > 5300000))
Platypus_ChIP_X5spec <- subset(Platypus_ChIP2, (Platypus_ChIP2$Chromosome == 'NC_041753.1' & Platypus_ChIP2$Start > 8500000 & Platypus_ChIP2$End < 69200000))

Platypus_ChIP_Xspec <- rbind(Platypus_ChIP_X1spec, Platypus_ChIP_X2spec, Platypus_ChIP_X3spec, Platypus_ChIP_X4spec, Platypus_ChIP_X5spec)

for(i in LISTOFHISTMODS){
  x <- sum(Platypus_ChIP_Xspec$Modification == i & Platypus_ChIP_Xspec$Sample == 'Female' & (Platypus_ChIP_Xspec$Chromosome == 'NC_041749.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041750.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041751.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041752.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041753.1'))
  write.csv(x, file = paste(i,'_female_Xspec_peaks.txt'))
  x <- sum(Platypus_ChIP_Xspec$Modification == i & Platypus_ChIP_Xspec$Sample == 'Male' & (Platypus_ChIP_Xspec$Chromosome == 'NC_041749.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041750.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041751.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041752.1' | Platypus_ChIP_Xspec$Chromosome == 'NC_041753.1'))
  write.csv(x, file = paste(i,'_male_Xspec_peaks.txt'))
}
