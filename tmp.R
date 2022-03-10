
rna <- read_delim("/g/data/xl04/pw8697/OAN/oan.normalised.counts.txt", delim = "\t")
tissues <- c("Fib","Hea","Liv")
prot <- list()
for (i in 1:length(tissues)){
    x <- read_delim(paste("/g/data/xl04/pw8697/Prot_ratios/OAN", tissues[i],".txt",sep=""), delim = "\t")
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
group_by(chr,ctype,type) %>% summarise(median = median(ratio, na.rm = T), x_start = min(start), x_end = max(start))

bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>%
select(chr,start,ratio=rFib) %>% mutate(type = "RNA"),
select(prot, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c("chrX1","chrX5")) %>%
mutate(ctype = case_when(chr == "chrX1" & start < 45000000 ~ "chrX1-PAR1", chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", TRUE ~ chr)) %>%
mutate(ratio = case_when(ratio > 4 ~ 4, TRUE ~ ratio)) %>%
ggplot(aes(x=start,y=log(ratio,2))) + geom_point(aes(color = type), size = 3) + 
geom_segment(data = medians, mapping = aes(x = x_start, xend = x_end, y = log(median,2), yend = log(median,2), color = type), size = 3) + 
facet_wrap(~chr,ncol=2, scales = "free_x") + theme_bw() + 
theme(text = element_text(size=20)) 

give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}


##for boxplots
autosomes <- paste("chr", c(1:21),sep="")

bind_rows(rna %>% mutate(rFib = fibroblastM/fibroblastF, rHea = heartM/heartF, rLiv = liverM / liverF) %>%
select(chr,geneid, start,ratio=rFib) %>% mutate(type = "RNA"),
select(prot, geneid, chr, start = pstart, ratio = pFib) %>% mutate(type = "Protein")) %>% filter(chr %in% c(autosomes, "chrX1", "chrX5")) %>%
mutate(ctype = case_when(chr == "chrX1" & start < 45000000 ~ "chrX1-PAR1",
                         chr == "chrX5" & start < 8500000 ~ "chrX5-PAR1", 
                         chr %in% autosomes ~ "autosome", TRUE ~ chr)) %>%
filter(!is.na(ratio)) %>%
group_by(chr,geneid) %>% mutate(count = n(), class = case_when(type == "RNA" & count == 2 ~ "Subset RNA", TRUE ~ type)) %>% 
ggplot(aes(x=ctype,y=log(ratio,2),fill = factor(class, levels = c("RNA","Subset RNA","Protein")))) + 
geom_boxplot(notch = T, outlier.shape = NA) +
stat_summary(fun.data = give.n, geom = "text", hjust=0.5, vjust=2, position = position_dodge(width = 0.75)) + 
theme_bw() + theme(text = element_text(size=20)) + coord_cartesian(ylim = c(-2, 2)) + xlim(c("autosome","chrX1-PAR1","chrX1","chrX5-PAR1","chrX5")) +
labs(fill="")
