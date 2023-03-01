##Set working directory
setwd("/Users/ashleymilton/Desktop/ChIP data local/Platypus_methylation")

#New male and female files

Platypus_RRBS_M <- read.delim("PlatypusRRBS_MaleTiles_50kb.txt", header = TRUE)
Platypus_RRBS_F <- read.delim("PlatypusRRBS_FemaleTiles_50kb.txt", header = TRUE)

Platypus_RRBS_M$percent_coverage <- ((Platypus_RRBS_M$numCs)/(Platypus_RRBS_M$coverage))*100
Platypus_RRBS_F$percent_coverage <- ((Platypus_RRBS_F$numCs)/(Platypus_RRBS_F$coverage))*100

Platypus_RRBS_M_small <- subset(Platypus_RRBS_M, select=-c(4, 5, 6, 7))
#Platypus_RRBS_M_small <- Platypus_RRBS_M_small[!(Platypus_RRBS_M_small$chr=="NC_000891.1"),]
Platypus_RRBS_F_small <- subset(Platypus_RRBS_F, select=-c(4, 5, 6, 7))

Platypus_RRBS_M_small$midpoint <- Platypus_RRBS_M_small$end - ((Platypus_RRBS_M_small$end - Platypus_RRBS_M_small$start)/2)
Platypus_RRBS_F_small$midpoint <- Platypus_RRBS_F_small$end - ((Platypus_RRBS_F_small$end - Platypus_RRBS_F_small$start)/2)


Platypus_RRBS_M_1 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041728.1"),]
Platypus_RRBS_M_2 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041729.1"),]
Platypus_RRBS_M_3 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041730.1"),]
Platypus_RRBS_M_4 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041731.1"),]
Platypus_RRBS_M_5 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041732.1"),]
Platypus_RRBS_M_6 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041733.1"),]
Platypus_RRBS_M_7 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041734.1"),]
Platypus_RRBS_M_8 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041735.1"),]
Platypus_RRBS_M_9 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041736.1"),]
Platypus_RRBS_M_10 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041737.1"),]
Platypus_RRBS_M_11 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041738.1"),]
Platypus_RRBS_M_12 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041739.1"),]
Platypus_RRBS_M_13 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041740.1"),]
Platypus_RRBS_M_14 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041741.1"),]
Platypus_RRBS_M_15 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041742.1"),]
Platypus_RRBS_M_16 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041743.1"),]
Platypus_RRBS_M_17 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041744.1"),]
Platypus_RRBS_M_18 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041745.1"),]
Platypus_RRBS_M_19 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041746.1"),]
Platypus_RRBS_M_20 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041747.1"),]
Platypus_RRBS_M_21 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041748.1"),]
Platypus_RRBS_M_X1 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041749.1"),]
Platypus_RRBS_M_X2 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041750.1"),]
Platypus_RRBS_M_X3 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041751.1"),]
Platypus_RRBS_M_X4 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041752.1"),]
Platypus_RRBS_M_X5 <- Platypus_RRBS_M_small[(Platypus_RRBS_M_small$chr=="NC_041753.1"),]


Platypus_RRBS_F_1 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041728.1"),]
Platypus_RRBS_F_2 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041729.1"),]
Platypus_RRBS_F_3 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041730.1"),]
Platypus_RRBS_F_4 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041731.1"),]
Platypus_RRBS_F_5 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041732.1"),]
Platypus_RRBS_F_6 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041733.1"),]
Platypus_RRBS_F_7 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041734.1"),]
Platypus_RRBS_F_8 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041735.1"),]
Platypus_RRBS_F_9 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041736.1"),]
Platypus_RRBS_F_10 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041737.1"),]
Platypus_RRBS_F_11 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041738.1"),]
Platypus_RRBS_F_12 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041739.1"),]
Platypus_RRBS_F_13 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041740.1"),]
Platypus_RRBS_F_14 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041741.1"),]
Platypus_RRBS_F_15 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041742.1"),]
Platypus_RRBS_F_16 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041743.1"),]
Platypus_RRBS_F_17 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041744.1"),]
Platypus_RRBS_F_18 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041745.1"),]
Platypus_RRBS_F_19 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041746.1"),]
Platypus_RRBS_F_20 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041747.1"),]
Platypus_RRBS_F_21 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041748.1"),]
Platypus_RRBS_F_X1 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041749.1"),]
Platypus_RRBS_F_X2 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041750.1"),]
Platypus_RRBS_F_X3 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041751.1"),]
Platypus_RRBS_F_X4 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041752.1"),]
Platypus_RRBS_F_X5 <- Platypus_RRBS_F_small[(Platypus_RRBS_F_small$chr=="NC_041753.1"),]


#Platypus_RRBS_M_small <- Platypus_RRBS_M_small[Platypus_RRBS_M_small$percent_coverage!=100 & Platypus_RRBS_M_small$percent_coverage!=0,]



library(ggplot2)
library(scales)


#Saving male/female together X1
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_X1.pdf", width=15, height=5)
ggplot()+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_X1$midpoint, y=Platypus_RRBS_F_X1$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_X1$midpoint, y=Platypus_RRBS_M_X1$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome X1", y = "% methylation coverage")+
  ggtitle("% methylation coverage for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#





#Saving male/female together X2
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_X2.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_X2$midpoint, y=GC_content_X2$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_X2$midpoint, y=GC_content_X2$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_X2$midpoint, y=Platypus_RRBS_F_X2$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_X2$midpoint, y=Platypus_RRBS_M_X2$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(0,100))+
  labs(x = "Position on Chromosome X2", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#




#Saving male/female together X3
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_X3.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_X3$midpoint, y=GC_content_X3$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_X3$midpoint, y=GC_content_X3$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_X3$midpoint, y=Platypus_RRBS_F_X3$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_X3$midpoint, y=Platypus_RRBS_M_X3$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(20,100))+
  labs(x = "Position on Chromosome X3", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#



#Saving male/female together X4
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_X4.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_X4$midpoint, y=GC_content_X4$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_X4$midpoint, y=GC_content_X4$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_X4$midpoint, y=Platypus_RRBS_F_X4$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_X4$midpoint, y=Platypus_RRBS_M_X4$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(0,100))+
  labs(x = "Position on Chromosome X4", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#



#Saving male/female together X5
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_X5.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_X5$midpoint, y=GC_content_X5$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_X5$midpoint, y=GC_content_X5$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_X5$midpoint, y=Platypus_RRBS_F_X5$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_X5$midpoint, y=Platypus_RRBS_M_X5$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(10,100))+
  labs(x = "Position on Chromosome X5", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#


#Saving male/female together 1
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_1.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_1$midpoint, y=GC_content_1$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_1$midpoint, y=GC_content_1$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_1$midpoint, y=Platypus_RRBS_F_1$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_1$midpoint, y=Platypus_RRBS_M_1$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 1", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#



#Saving male/female together 2
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_2.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_2$midpoint, y=GC_content_2$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_2$midpoint, y=GC_content_2$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_2$midpoint, y=Platypus_RRBS_F_2$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_2$midpoint, y=Platypus_RRBS_M_2$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 2", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 3
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_3.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_3$midpoint, y=GC_content_3$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_3$midpoint, y=GC_content_3$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_3$midpoint, y=Platypus_RRBS_F_3$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_3$midpoint, y=Platypus_RRBS_M_3$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 3", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 4
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_4.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_4$midpoint, y=GC_content_4$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_4$midpoint, y=GC_content_4$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_4$midpoint, y=Platypus_RRBS_F_4$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_4$midpoint, y=Platypus_RRBS_M_4$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 4", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 5
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_5.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_5$midpoint, y=GC_content_5$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_5$midpoint, y=GC_content_5$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_5$midpoint, y=Platypus_RRBS_F_5$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_5$midpoint, y=Platypus_RRBS_M_5$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 5", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#



#Saving male/female together 6
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_6.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_6$midpoint, y=GC_content_6$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_6$midpoint, y=GC_content_6$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_6$midpoint, y=Platypus_RRBS_F_6$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_6$midpoint, y=Platypus_RRBS_M_6$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 6", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 7
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_7.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_7$midpoint, y=GC_content_7$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_7$midpoint, y=GC_content_7$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_7$midpoint, y=Platypus_RRBS_F_7$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_7$midpoint, y=Platypus_RRBS_M_7$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 7", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 8
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_8.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_8$midpoint, y=GC_content_8$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_8$midpoint, y=GC_content_8$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_8$midpoint, y=Platypus_RRBS_F_8$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_8$midpoint, y=Platypus_RRBS_M_8$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 8", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 9
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_9.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_9$midpoint, y=GC_content_9$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_9$midpoint, y=GC_content_9$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_9$midpoint, y=Platypus_RRBS_F_9$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_9$midpoint, y=Platypus_RRBS_M_9$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 9", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 10
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_10.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_10$midpoint, y=GC_content_10$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_10$midpoint, y=GC_content_10$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_10$midpoint, y=Platypus_RRBS_F_10$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_10$midpoint, y=Platypus_RRBS_M_10$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 10", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#


#Saving male/female together 11
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_11.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_11$midpoint, y=GC_content_11$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_11$midpoint, y=GC_content_11$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_11$midpoint, y=Platypus_RRBS_F_11$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_11$midpoint, y=Platypus_RRBS_M_11$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 11", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 12
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_12.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_12$midpoint, y=GC_content_12$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_12$midpoint, y=GC_content_12$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_12$midpoint, y=Platypus_RRBS_F_12$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_12$midpoint, y=Platypus_RRBS_M_12$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 12", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 13
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_13.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_13$midpoint, y=GC_content_13$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_13$midpoint, y=GC_content_13$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_13$midpoint, y=Platypus_RRBS_F_13$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_13$midpoint, y=Platypus_RRBS_M_13$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 13", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 14
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_14.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_14$midpoint, y=GC_content_14$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_14$midpoint, y=GC_content_14$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_14$midpoint, y=Platypus_RRBS_F_14$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_14$midpoint, y=Platypus_RRBS_M_14$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 14", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 15
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_15.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_15$midpoint, y=GC_content_15$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_15$midpoint, y=GC_content_15$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_15$midpoint, y=Platypus_RRBS_F_15$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_15$midpoint, y=Platypus_RRBS_M_15$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 15", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 16
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_16.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_16$midpoint, y=GC_content_16$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_16$midpoint, y=GC_content_16$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_16$midpoint, y=Platypus_RRBS_F_16$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_16$midpoint, y=Platypus_RRBS_M_16$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 16", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#

#Saving male/female together 17
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_17.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_17$midpoint, y=GC_content_17$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_17$midpoint, y=GC_content_17$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_17$midpoint, y=Platypus_RRBS_F_17$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_17$midpoint, y=Platypus_RRBS_M_17$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 17", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#


#Saving male/female together 19
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_19.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_19$midpoint, y=GC_content_19$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_19$midpoint, y=GC_content_19$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_19$midpoint, y=Platypus_RRBS_F_19$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_19$midpoint, y=Platypus_RRBS_M_19$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 19", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()
#


#Saving male/female together 18
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_18.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_18$midpoint, y=GC_content_18$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_18$midpoint, y=GC_content_18$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_18$midpoint, y=Platypus_RRBS_F_18$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_18$midpoint, y=Platypus_RRBS_M_18$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 18", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()


#Saving male/female together 20
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_20.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_20$midpoint, y=GC_content_20$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_20$midpoint, y=GC_content_20$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_20$midpoint, y=Platypus_RRBS_F_20$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_20$midpoint, y=Platypus_RRBS_M_20$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 20", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()


#Saving male/female together 21 
pdf(file="/Users/ashleymilton/Library/CloudStorage/OneDrive-SharedLibraries-UNSW/Paul Waters - DC_Paper/Science/Paper_Figures/RRBS/Figures/Ash_methylation_030622/OAN_meth_GC_MF_21.pdf", width=15, height=5)
ggplot()+
  geom_line(aes(x=GC_content_21$midpoint, y=GC_content_21$gc*100, colour="darkorange1"))+
  geom_smooth(method = "loess", span = 0.03, size=0.5, aes(x=GC_content_21$midpoint, y=GC_content_21$gc*100, colour="gray30"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_F_21$midpoint, y=Platypus_RRBS_F_21$percent_coverage, colour="plum2", fill="plum2"))+
  geom_smooth(method = "loess", span = 0.05, aes(x=Platypus_RRBS_M_21$midpoint, y=Platypus_RRBS_M_21$percent_coverage, colour="cadetblue3", fill="cadetblue3"))+
  theme_classic()+
  scale_x_continuous(labels = comma)+
  coord_cartesian(ylim=c(25,100))+
  labs(x = "Position on Chromosome 21", y = "% methylation coverage")+
  ggtitle("% methylation coverage and GC content for platypus (50 kb tiles)")+
  scale_color_identity(name="Lines", labels=c("Male % methylation coverage LOBF", "GC content", "GC content LOBF", "Female % methylation coverage LOBF"))+
  scale_fill_identity(name="Fill", labels=c("Male % methylation coverage 95% CI", "Female % methylation coverage 95% CI"))
dev.off()

