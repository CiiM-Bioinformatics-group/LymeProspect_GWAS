library(ggplot2)
library(data.table)
library(tidyverse)
library(ggpubr)


#ORPlot

# Create labels
boxLabels = c("rs1061632_T Discovery", "rs1061632_T Validation")

rs2 <- fread(cmd = "egrep 'Z_STAT|ADD' out/lyme1_500FG.pheno.glm.logistic.hybrid | egrep 'Z_STAT|rs1061632'")
rs <- fread(cmd = "egrep 'Z_STAT|ADD' out/lyme2_300BCG.pheno.glm.logistic.hybrid | egrep 'Z_STAT|rs1061632'")

df <- data.frame(
  yAxis = length(boxLabels):1,
  boxOdds = c(1/rs$OR, 1/rs2$OR),
  boxCILow = c(1/rs$L95, 1/rs2$L95),
  boxCIHigh = c(1/rs$U95, 1/rs2$U95)
)

p <- ggplot(df, aes(x = boxOdds, y = yAxis))
p + geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), size = .5, height = .2, color = "gray50") +
  geom_point(size = 3.5, color = "orange") +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_y_continuous(breaks = df$yAxis, labels = boxLabels) +
  #scale_x_continuous(breaks = seq(1.0,3,3.1) ) +
  coord_trans(x = "log10") +
  ylab("") +
  xlab("Odds ratio (log scale)") + ggtitle("Susceptibility to Lyme disease")

ggsave("../SVG/Figure2.ORplot.svg", width = 8, height = 6)

#Mirror Locuszoomplot
chr = rs$`#CHROM`
min = rs$POS-5e5
max = rs$POS+5e5

#Get info
#Get info for both batches and calculate logP
locus1 <- fread(cmd = paste0("awk '{if($1 ==",chr,
                            " && $2 > ",min,
                            " && $2 < ",max,
                            ")print}' out/lyme1_500FG.pheno.glm.logistic.hybrid | grep ADD"))
colnames(locus1) <- colnames(rs)
locus1 <- locus1 %>% mutate(logP = log10(P))
locus1$batch <- "1"

#Add the imputation information
locus1_impu <- fread(cmd = paste0("awk -F ':' '{if($2 > ",min,
                                  " && $2 < ",max,
                                  ") print}' output/batch1/imputation/local/chr",chr, ".info"))
colnames(locus1_impu) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

locus1_impu <- locus1_impu %>% separate(SNPid, into = c('#CHROM', 'POS'))%>% dplyr::rename('REF' = Ref, 'ALT' = Alt)%>%
  mutate_at(c('#CHROM', 'POS'), as.integer)

locus1 <- left_join(locus1, locus1_impu, by = c('#CHROM', 'POS', 'REF', 'ALT'))

#Same for batch2, different sign of logP
locus2 <- fread(cmd = paste0("awk '{if($1 ==",chr,
                             " && $2 > ",min,
                             " && $2 < ",max,
                             ")print}' out/lyme2_300BCG.pheno.glm.logistic.hybrid | grep ADD"))
colnames(locus2) <- colnames(rs)
locus2 <- locus2 %>% mutate(logP = -log10(P))
locus2$batch <- "2"

#Add the imputation information
locus2_impu <- fread(cmd = paste0("awk -F ':' '{if($2 > ",min,
                                  " && $2 < ",max,
                                  ") print}' output/batch2/imputation/local/chr",chr, ".info"))
colnames(locus2_impu) <- c("SNPid","Ref","Alt","Alt_frq","MAF","AvgCall","Rsq","Genotyped","-","--","---","----","-----")

locus2_impu <- locus2_impu %>% separate(SNPid, into = c('#CHROM', 'POS'))%>% dplyr::rename('REF' = Ref, 'ALT' = Alt)%>%
  mutate_at(c('#CHROM', 'POS'), as.integer)

locus2 <- left_join(locus2, locus2_impu, by = c('#CHROM', 'POS', 'REF', 'ALT'))
#Merge both
locus <- rbind(locus1, locus2)
#Plotting
p1 <- ggplot(locus, aes(x=as.numeric(POS), y=logP))+
  geom_point(aes(color = Genotyped))+scale_color_manual(values = c('black', 'darkgrey'))+
  geom_text(data = locus[grepl('rs1061632',locus$ID),], aes(x=as.numeric(POS),y=logP+.5, label='rs1061632'), color = '#FAA51A', nudge_y = 0.2)+
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

eQTL <- c('ETV7', 'KCTD20', 'STK38', 'PI16', 'SRSF3')
#Plot genes
p3 <- ggplot(data = out.bm.genes.region) +
  geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, color = external_gene_name %in% eQTL)) +
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
ggsave("../SVG/Fig2.Locuszoom.svg", plot = plot, width = 7, height = 12)

#cQTL plots
#The effect of rs1061632 on other things
dosages <- as.data.frame(fread("zgrep -e rs1061632 -e HV 500FG/Genotype/merged/merged_dosages.txt.gz"))
dosages <- rename(dosages, 'V1' = id)
cytokines <- as.data.frame(fread("500FG/Molecular_phenotype/pheno_91cytokines_4Raul.csv"))
#A bit of having the same ppl
dosages <- dosages[,intersect(colnames(dosages), colnames(cytokines))]
cytokines <- cytokines[,intersect(colnames(dosages), colnames(cytokines))]


cQTL_boxplot <- function(rsID, Cyt, write = T){
  #Read genotype and cytokine
  print(Cyt)
  print(rsID)
  cQTL <- dosages[grepl(paste0(rsID,"$"), dosages$V1),]
  Cytokine_measured<- cytokines[cytokines$V1 == Cyt,]
  if(nrow(cQTL) > 1) stop("Too many snps")
  if(nrow(Cytokine_measured) > 1) stop("Too many snps")

  #Match column names
  Cytokine_measured <- Cytokine_measured[,colnames(cQTL)]

  #Bind
  df <- rbind(cQTL, Cytokine_measured)%>%
    remove_rownames()%>%
    column_to_rownames("V1")%>%
    t()%>%
    as.data.frame()

  colnames(df)[1] <-"snp"
  colnames(df)[2] <-"Cytokine"
  require(ggpubr)
  boxplot_theme <- function(){
    theme_bw()+theme(panel.grid = element_blank(), strip.background = element_blank(),axis.line = element_line(colour = "black"),
                     panel.border = element_blank(),panel.background = element_blank())
  }

  p <- ggboxplot(df, x = "snp",y = "Cytokine", color  = "snp",
                 add = c("point"), add.params = list(position = position_jitter(w = 0.05)))+
    ylab(paste("log2",Cyt))+xlab(rsID)+
    boxplot_theme()+scale_color_manual(values = c("deepskyblue4","deepskyblue", "lightskyblue3"))+
    stat_compare_means(comparisons = list(c(1,2),c(1,3),c(2,3)))
  if(write) ggsave(paste0("../SVG/Fig2.Boxplot_",Cyt,'_',rsID,".svg"), plot = p, width = 6, height = 6) else p
}

#First the cytokines
cQTL_boxplot('rs1061632', 'IL1b_Borreliamix_PBMC_24h')

#eQTLgen plots
#Read eQTLgen results for the snp
eQTLgen <- fread('zgrep -e rs1061632 -e SNPChr eqtlGen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
#Change to the assessed allele = T
eQTLgen <- mutate(eQTLgen, 'AssessedAllele' = 'T', 'OtherAllele' = 'C', Zscore = -Zscore)

ggplot(eQTLgen, aes(x=Zscore, y=reorder(GeneSymbol, Zscore)))+
  geom_col(aes(fill = as.factor(sign(Zscore))))+scale_fill_manual(values = c('darkred','darkblue'))+
  geom_text(aes(label = GeneSymbol, x = Zscore/2), color = 'white', fontface = "bold")+
  theme_minimal()+theme(panel.grid = element_blank(), axis.text.y = element_blank(), axis.line.x = element_line(), legend.position = 0)+
  ylab("Genes")+xlab("Z-score eQTLgen")
ggsave('../SVG/Fig2.eQTLgen.svg', width = 5, height = 5)

#RNAseq plots
pheno <- read.table("Raw/RNA/LB_phenotypes.txt", header = T)
counts <- read.table("Raw/RNA/LB_HV_normCounts_filt.txt", header = T)
colnames(counts) <- gsub("H.A.","H-A-", colnames(counts)) %>%
  gsub("LP.A.","LP-A-", .) %>%
  gsub("LP.R.", "LP-R-", .)

all(pheno[,1] == colnames(counts))

pheno_t <- pheno%>% mutate("sample" = gsub(".count.txt", "", sampleID)) %>%
  separate(remove = F, "sample", sep = "_", into = c("patient","s","t"))

#Gene symbols
library(clusterProfiler)
library(org.Hs.eg.db)
genesym <-  bitr(rownames(counts), fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = org.Hs.eg.db, drop = F)

#Get only genes I want
counts_sub_mat <- counts[genesym$ENSEMBL[genesym$SYMBOL %in%
                                       c('KCTD20','ETV7','PXT1', 'PI16','SRSF3','STK38','IRF7','RELN','CBS', 'IL10', 'IL1B', 'IL1RN', 'IL6', 'IFNG', 'IRF1', 'IRF3')],]
counts_sub <- counts_sub_mat %>% rownames_to_column('ENSEMBL')%>% as.data.frame() %>% reshape2::melt(variable.name = 'sampleID')
#Add gene symbol
genekey <- pull(genesym, SYMBOL, name = ENSEMBL)
counts_sub$gene <- recode(counts_sub$ENSEMBL, !!!genekey)

#1. Change in expression upon stimulation of KCTD20 and ETV7, and the eQTLS
eQTL <- c('ETV7', 'KCTD20', 'STK38', 'PI16', 'SRSF3')
counts_sub <- inner_join(counts_sub, pheno_t)%>%
  dplyr::select(-ENSEMBL)%>%
  pivot_wider(names_from = gene, values_from = value)

#Boxplots
for (i in eQTL){
  ggboxplot(counts_sub %>% subset(!grepl("t6",t) & s != 'S8')%>%mutate_at(i, log10), x="group.1", y = i, color= "s", palette = "jco", add = "jitter", line.color = "gray", add.params = list(position = position_jitter(w = 0.05)))+
    theme(axis.text.x = element_text(angle = 90))+scale_color_manual(values = c("deepskyblue4","lightskyblue3"))+
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       comparisons = list(c(1,2),c(3,4),c(1,3),c(2,4)))+
    ylab(paste('log10',i))
  ggsave(paste0("../SVG/Fig2.disease_healthy_t0_",i, ".svg"), width =7, height = 7)
}


#We find IRF7 has a TF motif that is distorted by the SNP
#Add IRF7 motif from https://regulomedb.org/regulome-search?regions=chr6%3A36457483-36457484&genome=GRCh37/thumbnail=motifs

#We know IRF7 is distorted, which genes is it affecting?
rownames(counts_sub_mat) <- recode(rownames(counts_sub_mat), !!!genekey)

#Clearly, ETV7 (or SRSF3)
svglite::svglite("../SVG/IRF7_correlation.svg", width = 6, height = 6)
ComplexHeatmap::Heatmap(cor(t(counts_sub_mat[c(eQTL, 'IRF7'),]), method = "spearman"), show_row_dend = F, show_column_dend = F)
dev.off()
#How significantly?
apply(counts_sub_mat[eQTL,], 1, function(x) cor.test(x,t(counts_sub_mat['IRF7',]), method ='spearman'))


#Calculating the ratio between both
counts_sub$ratio <- counts_sub$ETV7/counts_sub$IRF7
#How is rs1061632 affecting gene expression and the ratio between the two?
rs1061632 <- fread(cmd = 'grep -e LP -e rs1061632 input/lyme-merged_dosage.tsv')
rs1061632 <- column_to_rownames(rs1061632, 'V1')%>% t() %>% as.data.frame()%>%rownames_to_column('patient')
rs1061632$rs1061632 <- ifelse(rs1061632$`6:36457484:T:C%rs1061632` < 1.5, ifelse(rs1061632$`6:36457484:T:C%rs1061632` < 0.5, 0, 1), 2)
counts_sub <- left_join(counts_sub, rs1061632)
#Do the boxplots
#Boxplots
for (i in c(eQTL, 'ratio')){
  print(ggboxplot(counts_sub %>%mutate_at(i, log10) %>% na.omit(), x="rs1061632", color="group.1", y = i,  palette = "jco", add = "jitter", line.color = "gray", add.params = list(position = position_jitter(w = 0.05)))+
    facet_wrap(~group.1)+
    theme(axis.text.x = element_text(angle = 90))+
    stat_compare_means(method = "wilcox.test", label = "p.signif")+
    ylab(paste('log10',i)))
  #ggsave(paste0("../SVG/Fig2.Lyme_RNA_rs1061632_",i, ".svg"), width =5, height = 5)
}


#IRF7 hypothesis

#Reading the 500FG rnaseq data and finding the correlation between IRF7 and ETV7, checking if the correlation is different in healthy patients
irf7_gen <- fread("batch1/lyme1_500FG_rs1061632.traw")%>%
  column_to_rownames(var = "SNP")%>%
  dplyr::select(-c("CHR","(C)M","POS","COUNTED","ALT"))%>%
  t() %>% as.data.frame()%>% rownames_to_column(var = "patient")%>%
  separate(patient, into = 'patient')

#We are counting C, the alternative and minor allele, the major is T
irf7_gen <- mutate_at(irf7_gen, "rs1061632", function(x)
  ifelse(x == 0, "TT", ifelse(x == 1, "TC", ifelse(x == 2, "CC", x))))

irf7_counts <- fread("Raw/RNA/1508_Li_RNAseq.expression.genelevel.v75.htseq.txt")%>%
  column_to_rownames(var = "probe")

#Library size corrected counts
nnn <- DESeq2::DESeqDataSetFromMatrix(countData = irf7_counts, colData = data.frame("data" = colnames(irf7_counts), "yes" = "yes"), design = ~ 1)%>%
  DESeq2::estimateSizeFactors()

#Get the counts
irf7_rna_full <- DESeq2::counts(nnn, normalized = T)

#Get the genes
irf7_rna <- DESeq2::counts(nnn, normalized = T)%>%
  .[c("ENSG00000126456","ENSG00000010030","ENSG00000185507","ENSG00000126561",
      "ENSG00000112078","ENSG00000112079","ENSG00000112081","ENSG00000125347"),]%>%
  t() %>% as.data.frame()%>%rownames_to_column(var = "patient")%>%
  dplyr::rename("IRF7" = ENSG00000185507, "ETV7" = ENSG00000010030,
                "STAT5" = ENSG00000126561, "IRF3" = ENSG00000126456,
                "KCTD20" = ENSG00000112078, "STK38" = ENSG00000112079,
                "SRSF3" = ENSG00000112081, "IRF1" = ENSG00000125347)

#Merge genotype and RNA
irf7 <- left_join(irf7_rna, irf7_gen, "patient")

#Removing those three samples that have a 600 counts of ETV7
#irf7 <- irf7 %>% subset(ETV7 < 400)

#Correlation between two measures after removing one  outlier
ggplot(na.omit(irf7 %>% subset(ETV7 < 400)), aes(x = log10(IRF7), y = log10(ETV7), color = as.character(rs1061632)))+
  geom_point()+geom_smooth(method = "lm", se = T)+theme_minimal()
ggsave("../SVG/lm_irf7_500fg_snp_cut.svg")


lapply(c('ETV7','SRSF3','STK38', 'KCTD20','IRF7'), function(eGene){
  tmp <- irf7
  tmp$ratio <- tmp[,eGene]/tmp[,'IRF7']
  ggboxplot(na.omit(tmp), x= "rs1061632", y="ratio", color = "rs1061632", palette = "locuszoom", add = "jitter",
                  add.params = list(position = position_jitter(w = 0.05)))+scale_color_manual(values = c("deepskyblue4","deepskyblue", "lightskyblue3"))+
    stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2)))+
    ylab(paste0('ratio ', eGene,'/IRF7'))+theme(legend.position = 0)
  ggsave(paste0("../SVG/Figure3.irf7_",eGene,"_ratio_boxplot.svg"), width = 6, height = 6)

  tmp$logexpr <- log10(tmp[,eGene])
  ggboxplot(na.omit(tmp), x= "rs1061632", y='logexpr', color = "rs1061632", palette = "locuszoom", add = "jitter", add.params = list(position = position_jitter(w = 0.05)))+scale_color_manual(values = c("deepskyblue4","deepskyblue", "lightskyblue3"))+
    stat_compare_means(method = "wilcox.test", comparisons = list(c(1,2)))+
    ylab(paste0('log10 ', eGene))+theme(legend.position = 0)
  ggsave(paste0("../SVG/Figure3.",eGene,"_expr_boxplot.svg"), width = 6, height = 6)

})
