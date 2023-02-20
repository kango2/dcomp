library(tidyverse)

rna <- read_delim("~/Downloads/dcomp-main/data/OAN/oan.normalised.counts.txt", delim = "\t")
tissues <- c("Fib","Hea","Liv")
prot <- list()
for (i in 1:length(tissues)){
    x <- read_delim(paste("~/Downloads/dcomp-main/data/Prot_ratios/OAN", tissues[i],".txt",sep=""), delim = "\t")
    x <- x %>% mutate(chromosome_name = paste("chr",chromosome_name,sep=""), tissue = tissues[i])
    colnames(x) <- c("chr", "pstart","pend","pstrand", "refseqacc", "uniprotid", "geneid","gmeanratio", "tissue")
    prot[[i]] <- x
}
prot <- bind_rows(prot)
prot <- prot %>% pivot_wider(names_from = tissue, values_from = gmeanratio)
colnames(prot)[8:10] <- c("pFib", "pHea", "pLiv")

medians <- bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>%
select(chr,start,ratio=rFib) %>% mutate(type = "RNA"),
select(prot, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chrX1","chrX5")) %>%
mutate(ctype = case_when(chr == "chrX1" & start < 45000000 ~ "chrX1-PAR1", chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", TRUE ~ chr)) %>%
  mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
group_by(chr,ctype,type) %>% summarise(median = median(ratio, na.rm = T), x_start = min(start), x_end = max(start)) %>% 
  mutate(colour=case_when(type=="RNA" ~ "RNAmed")) %>%
  mutate(colour=case_when(type=="Protein" ~ "Protmed", TRUE ~ colour))

#chromosome plot
colours <- c("dodgerblue", "dodgerblue3","red3","red")

OANFibChrplot1 <- bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>%
select(chr,start,ratio=rFib) %>% mutate(type = "RNA"),
select(prot, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chrX1","chrX5")) %>%
mutate(ctype = case_when(chr == "chrX1" & start < 45000000 ~ "chrX1-PAR1", chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", TRUE ~ chr)) %>%
mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
ggplot(aes(x=start,y=log(ratio,2))) + geom_point(aes(color = type), size = 0.75, alpha = 0.5, shape = 18) + 
geom_segment(data = medians, mapping = aes(x = x_start, xend = x_end, y = log(median,2), yend = log(median,2), color = colour), size = 0.75,) + 
facet_wrap(~chr,ncol=2, scales = "free_x") + theme_bw() + 
theme(text = element_text(size=20)) +
coord_cartesian(ylim = c(-2, 2)) +
scale_color_manual(values = colours)

pdf("~/Downloads/OANFibChrplot1.pdf", height=2.5, width=11)
OANFibChrplot1
dev.off()

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

##for boxplots
autosomes <- paste("chr", c(1:21),sep="")

OANFibplot1 <- bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>% 
            select(chr,geneid, start,ratio=rFib) %>% mutate(type = "RNA"), 
          select(prot, geneid, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c(autosomes, "chrX1", "chrX5")) %>% 
  mutate(ctype = case_when(chr == "chrX1" & start < 45000000 ~ "chrX1-PAR1",
                         chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", 
                         chr %in% autosomes ~ "autosome", TRUE ~ chr)) %>%
filter(!is.na(ratio)) %>%
group_by(chr,geneid) %>% mutate(count = n(), class = case_when(type == "RNA" & count == 2 ~ "Subset RNA", TRUE ~ type)) %>% 
ggplot(aes(x=ctype,y=log(ratio,2),fill = factor(class, levels = c("RNA","Subset RNA","Protein")))) + 
geom_boxplot(notch = T, outlier.shape = NA) +
#stat_summary(fun.data = give.n, geom = "text", hjust=0.5, vjust=2, position = position_dodge(width = 0.75)) +
theme_bw() + theme(text = element_text(size=20)) + coord_cartesian(ylim = c(-2.5, 2)) + xlim(c("autosome","chrX1-PAR1","chrX1","chrX5-PAR1","chrX5")) +
labs(fill="")

pdf("~/Downloads/OANFibplot1.pdf", height=4, width=12)
OANFibplot1
dev.off()

library(rcompanion)

#calculate p values
ratios <- bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>%
                       select(chr,geneid,start,ratio=rFib) %>% mutate(type = "RNA"),
                     select(prot, geneid, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c(autosomes,"chrX1","chrX5")) %>%
  mutate(ctype = case_when(chr == "chrX1" & start < 45000000 ~ "chrX1-PAR1", chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", chr %in% autosomes ~ "autosome", TRUE ~ chr)) %>%
  filter(!is.na(ratio)) %>%
  group_by(chr,geneid) %>% mutate(count = n(), class = case_when(type == "RNA" & count == 2 ~ "Subset RNA", TRUE ~ type)) %>% 
  mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
  group_by(ctype,type,class) %>% 
  unite(ctype, ctype, class, sep = "_", remove=FALSE)


p_adjust_values_withautosomes <- pairwiseMedianTest(ratio ~ ctype,
                   data   = ratios,
                   method = "fdr")

ratios_RNA <- filter(ratios, class== "RNA")
p_adjust_RNA <- pairwiseMedianTest(ratio ~ ctype,
                                      data   = ratios_RNA,
                                      method = "fdr")

ratios_subsetRNA <- filter(ratios, class== "Subset RNA")
p_adjust_subsetRNA <- pairwiseMedianTest(ratio ~ ctype,
                                   data   = ratios_subsetRNA,
                                   method = "fdr")

ratios_Protein <- filter(ratios, class== "Protein")
p_adjust_Protein <- pairwiseMedianTest(ratio ~ ctype,
                                   data   = ratios_Protein,
                                   method = "fdr")


#Make p values data into excel spreadsheet
library("writexl")
write_xlsx(p_adjust_values_withautosomes, "~/Downloads/Platypus_Fib_p_values_with_autosomes.xlsx")



OANFibplotctype <- bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>% 
            select(chr,geneid, start,ratio=rFib) %>% mutate(type = "RNA"), 
          select(prot, geneid, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chrX1", "chrX2", "chrX3", "chrX4", "chrX5")) %>% 
  mutate(ctype = case_when(chr == "chrX1" & start < 45100000 ~ "PAR",
                           chr == "chrX2" & start < 5700000  ~ "PAR",
                           chr == "chrX2" & start > 20000000  ~ "PAR",
                           chr == "chrX3" & start < 18500000  ~ "PAR",
                           chr == "chrX3" & start > 31900000  ~ "PAR",
                           chr == "chrX4" & start < 5300000   ~ "PAR",
                           chr == "chrX5" & start < 8500000 ~ "PAR", 
                           chr == "chrX5" & start > 69200000 ~ "PAR", 
                           TRUE ~ chr)) %>% 
  mutate(ctype = case_when(grepl("chr", ctype) ~ "X", TRUE ~ ctype)) %>%
  filter(!is.na(ratio)) %>%
  group_by(chr,geneid) %>% mutate(count = n(), class = case_when(type == "RNA" & count == 2 ~ "Subset RNA", TRUE ~ type)) %>% 
  ggplot(aes(x=ctype,y=log(ratio,2),fill = factor(class, levels = c("RNA","Subset RNA","Protein")))) + 
  geom_boxplot(notch = T, outlier.shape = NA) +
  stat_summary(fun.data = give.n, geom = "text", hjust=0.5, vjust=2, position = position_dodge(width = 0.75)) +
  theme_bw() + theme(text = element_text(size=20)) + coord_cartesian(ylim = c(-2.5, 2)) + xlim(c("PAR","X")) +
  labs(fill="")

pdf("~/Downloads/OANFibplotctype.pdf", height=4, width=6)
OANFibplotctype
dev.off()

#to get median numbers on plot
#give.n <- function(x){
  #return(c(y = median(x)*1.05, label = median(x))) 
  # experiment with the multiplier to find the perfect position
#}


ratios <- bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>%
                      select(chr,geneid,start,ratio=rFib) %>% mutate(type = "RNA"),
                    select(prot, geneid, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chrX1", "chrX2", "chrX3", "chrX4", "chrX5")) %>%
  mutate(ctype = case_when(chr == "chrX1" & start < 45100000 ~ "PAR",
                           chr == "chrX2" & start < 5700000  ~ "PAR",
                           chr == "chrX2" & start > 20000000  ~ "PAR",
                           chr == "chrX3" & start < 18500000  ~ "PAR",
                           chr == "chrX3" & start > 31900000  ~ "PAR",
                           chr == "chrX4" & start < 5300000   ~ "PAR",
                           chr == "chrX5" & start < 8500000 ~ "PAR", 
                           chr == "chrX5" & start > 69200000 ~ "PAR", 
                           TRUE ~ chr)) %>%
  mutate(ctype = case_when(grepl("chr", ctype) ~ "X", TRUE ~ ctype)) %>%
  filter(!is.na(ratio)) %>%
  group_by(chr,geneid) %>% mutate(count = n(), class = case_when(type == "RNA" & count == 2 ~ "Subset RNA", TRUE ~ type)) %>% 
  mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
  group_by(ctype,type,class) %>% 
  unite(ctype, ctype, class, sep = "_", remove=FALSE)


p_adjust_values_withautosomes <- pairwiseMedianTest(ratio ~ ctype,
                                                    data   = ratios,
                                                    method = "fdr")


chr == "chrX1" & start < 45100000 ~ "chrX1-PAR1",
chr == "chrX2" & start < 5700000  ~ "chrX2-PAR1",
chr == "chrX2" & start > 20000000  ~ "chrX2-PAR2",
chr == "chrX3" & start < 18500000  ~ "chrX3-PAR1",
chr == "chrX3" & start > 31900000  ~ "chrX3-PAR2",
chr == "chrX4" & start < 5300000   ~ "chrX4-PAR1",
chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", 
chr == "chrX5" & start > 69200000 ~ "chrX5-PAR2"
