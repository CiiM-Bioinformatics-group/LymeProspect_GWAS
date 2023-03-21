library(ggplot2)
library(data.table)
library(tidyverse)
library(ggpubr)


rs <- fread(cmd = "egrep 'Z_STAT|ADD' /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_500FG.pheno.glm.logistic.hybrid | grep -w -e Z_STAT -e rs4110197")

#Mirror Locuszoomplot
chr = rs$`#CHROM`
min = rs$POS-2.5e5
max = rs$POS+2.5e5

#Get info
#Get info for both batches and calculate logP
locus1 <- fread(cmd = paste0("awk '{if($1 ==",chr,
                             " && $2 > ",min,
                             " && $2 < ",max,
                             ")print}' /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_500FG.pheno.glm.logistic.hybrid | grep ADD"))
colnames(locus1) <- colnames(rs)
locus1 <- locus1 %>% mutate(logP = -log10(P))
locus1$batch <- "1"

#Add the imputation information
locus1_impu <- fread(cmd = paste0("awk -F ':' '{if($2 > ",min,
                                  " && $2 < ",max,
                                  ") print}' /vol/projects/CIIM/Lyme_GWAS/GWAS/Genetics/output/batch1/imputation/local/chr",chr, ".info"))
colnames(locus1_impu) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

locus1_impu <- locus1_impu %>% separate(SNPid, into = c('#CHROM', 'POS'))%>% dplyr::rename('REF' = Ref, 'ALT' = Alt)%>%
  mutate_at(c('#CHROM', 'POS'), as.integer)

locus1 <- left_join(locus1, locus1_impu, by = c('#CHROM', 'POS', 'REF', 'ALT'))

#Same for batch2, different sign of logP
locus2 <- fread(cmd = paste0("awk '{if($1 ==",chr,
                             " && $2 > ",min,
                             " && $2 < ",max,
                             ")print}' /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme2_300BCG.pheno.glm.logistic.hybrid | grep ADD"))
colnames(locus2) <- colnames(rs)
locus2 <- locus2 %>% mutate(logP = log10(P))
locus2$batch <- "2"

#Add the imputation information
locus2_impu <- fread(cmd = paste0("awk -F ':' '{if($2 > ",min,
                                  " && $2 < ",max,
                                  ") print}' /vol/projects/CIIM/Lyme_GWAS/GWAS/Genetics/output/batch2/imputation/local/chr",chr, ".info"))
colnames(locus2_impu) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

locus2_impu <- locus2_impu %>% separate(SNPid, into = c('#CHROM', 'POS'))%>% dplyr::rename('REF' = Ref, 'ALT' = Alt)%>%
  mutate_at(c('#CHROM', 'POS'), as.integer)

locus2 <- left_join(locus2, locus2_impu, by = c('#CHROM', 'POS', 'REF', 'ALT'))
#Merge both
locus <- rbind(locus1, locus2)
locus <- locus %>% separate(ID, into =c('#CHROM','POS','REF','ALT','rsID'), remove = F)
#Plotting
p1 <- ggplot(locus, aes(x=as.numeric(POS), y=logP))+
  geom_point(aes(color = Genotyped))+scale_color_manual(values = c('black', 'darkgrey'))+
    ggrepel::geom_text_repel(data = locus[grepl('rs4110197|rs2232950',locus$ID),], aes(x=as.numeric(POS),y=logP, label=rsID), color = '#FAA51A', min.segment.length = 0,
                           box.padding = 0.5)+
  geom_hline(aes(yintercept = 0))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab(paste0("-log10(Pvalue)"))+xlab("")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Adding genes
library(biomaRt)
#Snp reference
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) # we will need an additional mart for genes
#extract genes
out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'),
  filters = c('chromosome_name','start','end'),
  values = list(chr, min, max),
  mart = gene.ensembl)
## define plot range for x-axis
plot.range <- c(min(min, out.bm.genes.region$start_position),
                max(max, out.bm.genes.region$end_position))

## rank gene_biotype label
out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype),
                                                                                     "protein_coding"), external_gene_name = fct_reorder2(external_gene_name,
                                                                                                                                          start_position, gene_biotype_fac, .desc = TRUE))%>%
  filter(gene_biotype == "protein_coding")

#Plot genes
p3 <- ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, color = 'SCGB1D2')) +
  scale_color_manual(values = c('black','#FAA51A'))+
  coord_flip() + ylab("") +
  ylim(plot.range) +
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name),
            fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) +
  #labs(caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
  theme_minimal()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.text.x = element_text(angle = 90),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position=0,
        panel.grid = element_blank(), axis.line.x = element_line()) +
  expand_limits(y=c(-1, 1))
library(patchwork)
plot <- p1/p3+plot_layout(heights = c(3,1))
plot
ggsave("../SVG/Fig6.Locuszoom_mirror.svg", plot = plot, width = 10, height = 10)


#Non-mirror locuszoom

#Get info
#Get info for both batches and calculate logP
locus <- fread(cmd = paste0("awk -F: '{if($1 ==",chr,
                             " && $2 > ",min,
                             " && $2 < ",max,
                             ")print}' /vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_lyme2.metal.tbl"))
colnames(locus) <- c('ID','REF','ALT','beta','se','P', 'direction')
locus <- locus %>% mutate(logP = -log10(P))

#Add the imputation information
locus_impu <- fread(cmd = paste0("awk -F ':' '{if($2 > ",min,
                                  " && $2 < ",max,
                                  ") print}' /vol/projects/CIIM/Lyme_GWAS/GWAS/Genetics/output/batch1/imputation/local/chr",chr, ".info"))
colnames(locus_impu) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

locus_impu <- locus_impu %>% separate(SNPid, into = c('#CHROM', 'POS'))%>% dplyr::rename('REF' = Ref, 'ALT' = Alt)

locus <- locus %>% separate(ID, into =c('#CHROM','POS','REF','ALT','rsID'), remove = F)

locus <- left_join(locus, locus_impu, by = c('#CHROM', 'POS', 'REF', 'ALT'))

#Plotting
p1 <- ggplot(locus, aes(x=as.numeric(POS), y=logP))+
  geom_point(aes(color = Genotyped))+scale_color_manual(values = c('black', 'darkgrey'))+
  ggrepel::geom_text_repel(data = locus[grepl('rs4110197|rs2232950',locus$ID),], aes(x=as.numeric(POS),y=logP, label=rsID), color = '#FAA51A', min.segment.length = 0,
                           box.padding = 0.5)+
  geom_hline(aes(yintercept = 0))+
  theme_bw()+theme(panel.grid = element_blank())+
  ylab(paste0("-log10(Pvalue)"))+xlab("")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#Adding genes
library(biomaRt)
#Snp reference
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37) # we will need an additional mart for genes
#extract genes
out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'),
  filters = c('chromosome_name','start','end'),
  values = list(chr, min, max),
  mart = gene.ensembl)
## define plot range for x-axis
plot.range <- c(min(min, out.bm.genes.region$start_position),
                max(max, out.bm.genes.region$end_position))

## rank gene_biotype label
out.bm.genes.region <- out.bm.genes.region %>% mutate(gene_biotype_fac = fct_relevel(as.factor(gene_biotype),
                                                                                     "protein_coding"), external_gene_name = fct_reorder2(external_gene_name,
                                                                                                                                          start_position, gene_biotype_fac, .desc = TRUE))%>%
  filter(gene_biotype == "protein_coding")

#Plot genes
p3 <- ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, color = 'SCGB1D2')) +
  scale_color_manual(values = c('black','#FAA51A'))+
  coord_flip() + ylab("") +
  ylim(plot.range) +
  geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name),
            fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) +
  #labs(caption = paste0("Data source: ", gene.ensembl@host, " + Data set: ", gene.ensembl@dataset), color = "Gene Biotype") +
  theme_minimal()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),axis.text.x = element_text(angle = 90),
        axis.ticks.y=element_blank(),
        strip.text.y = element_text(angle = 0),
        legend.position=0,
        panel.grid = element_blank(), axis.line.x = element_line()) +
  expand_limits(y=c(-1, 1))
library(patchwork)
plot <- p1/p3+plot_layout(heights = c(3,1))
plot
ggsave("../SVG/Fig6.Locuszoom_meta_analysis.svg", plot = plot, width = 10, height = 10)
