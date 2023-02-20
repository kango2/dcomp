library(tidyverse)

rna <- read_delim("~/Downloads/dcomp-main/data/OAN/mmu.normalised.counts.txt", delim = "\t")
tissues <- c("Hea","Liv")
prot <- list()
for (i in 1:length(tissues)){
    x <- read_delim(paste("~/Downloads/dcomp-main/data/Prot_ratios/MMU", tissues[i],".txt",sep=""), delim = "\t")
    x <- x %>% mutate(chromosome_name = paste("chr",chromosome_name,sep=""), tissue = tissues[i])
    colnames(x) <- c("chr", "pstart","pend","pstrand", "refseqacc", "uniprotid", "geneid","gmeanratio", "tissue")
    prot[[i]] <- x
}
prot <- bind_rows(prot)
prot <- prot %>% pivot_wider(names_from = tissue, values_from = gmeanratio)
colnames(prot)[8:9] <- c("pHea", "pLiv")

colnames(rna)[1] <- c("geneid")

medians <- bind_rows(rna %>% mutate(rHea = ((mhm1+mhm2)/2)/((mhf1+mhf2)/2), rLiv = ((mlm1+mlm2)/2)/((mlf1+mlf2)/2)) %>%
select(chr,start=gstart,ratio=rHea) %>% mutate(type = "RNA"),
select(prot, chr, start=pstart, ratio = pHea) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chr3","chrX")) %>%
mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
group_by(chr,type) %>% summarise(median = median(ratio, na.rm = T), x_start = min(start), x_end = max(start)) %>% 
  mutate(colour=case_when(type=="RNA" ~ "RNAmed")) %>%
  mutate(colour=case_when(type=="Protein" ~ "Protmed", TRUE ~ colour))

#chromosome plot
colours <- c("dodgerblue", "dodgerblue3","red3","red")

MMUHeaChrplot1 <- bind_rows(rna %>% mutate(rHea = ((mhm1+mhm2)/2)/((mhf1+mhf2)/2), rLiv = ((mlm1+mlm2)/2)/((mlf1+mlf2)/2)) %>%
select(chr,start=gstart,ratio=rHea) %>% mutate(type = "RNA"),
select(prot, chr, start=pstart, ratio = pHea) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chr3","chrX")) %>%
mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
ggplot(aes(x=start,y=log(ratio,2))) + geom_point(aes(color = type), size = 0.75, alpha = 0.5, shape = 18) + 
geom_segment(data = medians, mapping = aes(x = x_start, xend = x_end, y = log(median,2), yend = log(median,2), color = colour), size = 0.75) + 
facet_wrap(~chr,ncol=2, scales = "free_x") + theme_bw() + 
theme(text = element_text(size=20)) +
coord_cartesian(ylim = c(-2, 2)) +
scale_color_manual(values = colours) 

pdf("~/Downloads/MMUHeaChrplot1.pdf", height=2.5, width=11)
MMUHeaChrplot1
dev.off()

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


##for boxplots
autosomes <- paste("chr", c(1:19),sep="")

MMUHeaplot1 <- bind_rows(rna %>% mutate(rHea = ((mhm1+mhm2)/2)/((mhf1+mhf2)/2), rLiv = ((mlm1+mlm2)/2)/((mlf1+mlf2)/2)) %>% 
            select(chr,geneid, start=gstart,ratio=rHea) %>% mutate(type = "RNA"), 
          select(prot, geneid, chr, start=pstart, ratio=pHea) %>% mutate(type = "Protein")) %>% filter(chr %in% c(autosomes, "chrX")) %>% 
  mutate(ctype = case_when(chr %in% autosomes ~ "autosome", TRUE ~ chr)) %>%
filter(!is.na(ratio)) %>%
group_by(chr,geneid) %>% mutate(count = n(), class = case_when(type == "RNA" & count == 2 ~ "Subset RNA", TRUE ~ type)) %>% 
ggplot(aes(x=ctype,y=log(ratio,2),fill = factor(class, levels = c("RNA","Subset RNA","Protein")))) + 
geom_boxplot(notch = T, outlier.shape = NA) +
#stat_summary(fun.data = give.n, geom = "text", hjust=0.5, vjust=2, position = position_dodge(width = 0.75)) + 
theme_bw() + theme(text = element_text(size=20)) + coord_cartesian(ylim = c(-2.5, 2.5)) + xlim(c("autosome","chrX")) +
labs(fill="")
pdf("~/Downloads/MMUHeaplot1.pdf", height=4, width=12)
MMUHeaplot1
dev.off()

library(rcompanion)

#calculate p values
ratios <- bind_rows(rna %>% mutate(rHea = ((mhm1+mhm2)/2)/((mhf1+mhf2)/2), rLiv = ((mlm1+mlm2)/2)/((mlf1+mlf2)/2)) %>%
                      select(chr,geneid,start=gstart,ratio=rHea) %>% mutate(type = "RNA"),
                    select(prot, geneid, chr, start = pstart, ratio = pHea) %>% mutate(type = "Protein")) %>% filter(chr %in% c(autosomes,"chrX")) %>%
  mutate(ctype = case_when(chr %in% autosomes ~ "autosome", TRUE ~ chr)) %>%
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
write_xlsx(p_adjust_values_withautosomes, "~/Downloads/MMU_Hea_p_values_with_autosomes.xlsx")
