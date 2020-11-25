library(tidyverse)
library(edgeR)

setwd("/media/dstore/RNASeq_DC_Data/")

species <- c("mmu", "mdo", "oan", "gga")

##chromosome names map
datalist = list()
for (i in 1:length(species)) {
  file <- list.files(paste("refdata/", species[i], sep = ""), "*_assembly_report.txt", full.names = T)
  tmp <- read_delim(file, "\t", comment = "#", col_names = F) %>% 
    dplyr::select(chr = X1, acc = X7, len = X9) %>% 
    mutate(chr = str_remove(chr, "Super_Scaffold_")) %>%
    mutate(species = species[i], chr = paste("chr", chr, sep = ""))
  datalist[[i]] <- tmp
}
chrnamemap <- do.call(bind_rows, datalist)

##sample level metadata table

metadata <- read_delim("metadata.txt", "\t") %>% mutate(species = tolower(species)) 

##summary information about samples
#0 (unstranded), 1 (stranded) and 2 (reversely stranded)
#gga and mmu = s2 strand orientation for RNA
#mdo = s1 strand orientation for RNA
#oan = s0 strand orientation for RNA

files <- list.files(".", "*.txt.summary", full.names = T)
datalist <- list()
for ( i in 1:length(files)){
  tmp <- read_delim(files[i],"\t") %>% pivot_longer(!Status, names_to = "filenames", values_to = "counts")
  datalist[[i]] <- tmp
}

countsummary <- do.call(bind_rows, datalist)
countsummary <- separate(countsummary, col = filenames, into = c("bampath", "sampleid"), sep = ":")
countsummary <- pivot_wider(countsummary, names_from = "Status", values_from = "counts")
##merge metadata with countsummary
metadata <- left_join(metadata, countsummary, by = "sampleid")
write_delim(metadata, "metadata.counts.txt", delim = "\t")

##library size for normalisation later
libsize <- dplyr::group_by(metadata,scode) %>% dplyr::summarise(libsize=sum(readcounts))

processcounts <- function(counts = counts, species = "gga"){
  geneinfo <- dplyr::select(counts, Geneid, Chr, Start, End, Strand, Length)
  counts <- dplyr::select(counts, -Chr, -Start, -End, -Strand, -Length)
  colnames(counts) <- c("geneid", str_split(colnames(counts), ":", simplify = T)[,2][-1])
  counts <- counts %>% pivot_longer(!geneid, names_to = "sampleid", values_to = "counts")
  counts <- left_join(counts, metadata, by="sampleid") %>% 
    dplyr::select(geneid, sampleid, scode, counts) %>%
    dplyr::group_by(geneid,scode) %>% 
    dplyr::summarise(counts = sum(counts)) %>%
    pivot_wider(names_from = scode, values_from = counts)

  write_delim(counts, paste(species,".composite.counts.txt",sep = ""), delim = "\t")

  geneinfo <- left_join(mutate(geneinfo, acc = str_split(Chr,";",simplify = T)[,1],
                               gstart = str_split(Start,";",simplify = T)[,1],
                               gend = str_match(geneinfo$End,"\\S+;(\\S+)")[,2],
                               gstrand = str_split(Strand,";",simplify = T)[,1]),
                        chrnamemap, by = "acc")
  write_delim(geneinfo, paste(species,".geneinfo.txt",sep = ""), delim = "\t")

  y <- DGEList(counts=counts[,-1], 
              group = str_extract(colnames(counts[,-1]), "\\S{3}"), 
              lib.size = filter(libsize, scode %in% colnames(counts)) %>% pull(libsize),
              genes = left_join(dplyr::select(counts, Geneid=geneid), geneinfo, by = "Geneid"),
              remove.zeros = T
              )
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=F]
  y <- calcNormFactors(y)

  ##normalised counts
  ncounts <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(cpm(y)))
  lncounts <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(cpm(y, log = T)))
  rpkm <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(rpkm(y)))
  lrpkm <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(rpkm(y, log = T)))
  gncounts <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(cpmByGroup(y)))
  glncounts <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(cpmByGroup(y, log = T)))
  grpkm <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(rpkmByGroup(y)))
  glrpkm <- bind_cols(y$genes[,c("Geneid","species","chr","gstart","gend","gstrand","Length")], data.frame(rpkmByGroup(y, log = T)))

  write_delim(ncounts, paste(species,".normalised.counts.txt",sep = ""), delim = "\t")
  write_delim(lncounts, paste(species,".normalised.logcounts.txt",sep = ""), delim = "\t")
  write_delim(rpkm, paste(species,".rpkm.counts.txt",sep = ""), delim = "\t")
  write_delim(lrpkm, paste(species,".rpkm.logcounts.txt",sep = ""), delim = "\t")
  write_delim(gncounts, paste(species,".gnormalised.counts.txt",sep = ""), delim = "\t")
  write_delim(glncounts, paste(species,".gnormalised.logcounts.txt",sep = ""), delim = "\t")
  write_delim(grpkm, paste(species,".grpkm.counts.txt",sep = ""), delim = "\t")
  write_delim(glrpkm, paste(species,".grpkm.logcounts.txt",sep = ""), delim = "\t")
}

##read in counts for each species
##mean count threshold
##gga counts
counts <- read_delim("ggacounts.s2.txt","\t",col_names = T, comment = "#")
processcounts(counts=counts, species="gga")

counts <- read_delim("mmucounts.s2.txt", "\t", col_names = T, comment = "#")
processcounts(counts=counts, species="mmu")

counts <- read_delim("mdocounts.s1.txt", "\t", col_names = T, comment = "#")
processcounts(counts=counts, species="mdo")

counts <- read_delim("oancounts.fibroblast.s0.txt", "\t", col_names = T, comment = "#")
counts2 <- read_delim("oancounts.s0.txt", "\t", col_names = T, comment = "#")
counts <- left_join(counts, dplyr::select(counts2,-Chr,-Start,-End,-Strand,-Length), by = "Geneid")
processcounts(counts=counts, species="oan")

