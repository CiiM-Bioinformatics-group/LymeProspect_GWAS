#!/usr/bin/env Rscript
require("ggman")
require("ggplot2")
require("tidyverse")
require("data.table")

.convert2posX <- function(chr, bp, build) {
  if (length(chr) != length(bp)) {
    stop("SIZE DIFFER.");
  }
  
  n_chr = length(unique(chr))
  if (n_chr == 1) {
    return (list(posX = bp,
                 breaks = ggplot2::waiver(),
                 labels = ggplot2::waiver(),
                 xlabel = paste("Chromosome", chr[1])))
  }
  
  max_chr = max(chr)
  if (build == "hg19") {
    return (list(posX = ..convert2posX(chr, bp, .hg19),
                 breaks = .hg19.breaks[1:max_chr],
                 labels = .labels[1:max_chr],
                 xlabel = "Chromosome"))
  } else {
    posX = numeric(max_chr)
    breaks = numeric(max_chr)
    for (i in 1:max_chr) {
      size = max(bp[chr == i]) * 1.1
      posX[i + 1] = posX[i] + size
      breaks[i] = posX[i] + size / 2
    }
    return (list(posX = ..convert2posX(chr, bp, posX),
                 breaks = breaks[1:max_chr],
                 labels = .labels[1:max_chr],
                 xlabel = "Chromosome"))
  }
}

..convert2posX <- function (chr, bp, posX){
  .Call("ggman_convert2posX", PACKAGE = "ggman", chr, bp, posX)
}
.labels = c(as.character(seq(22)), "X", "Y", "XY", "MT")

# hg19 chromosome sizes
.hg19 = c(0,           # chr1
          249250621,
          492449994,
          690472424,
          881626700,
          1062541960,
          1233657027,
          1392795690,
          1539159712,
          1680373143,  # chr10
          1815907890,
          1950914406,
          2084766301,
          2199936179,
          2307285719,
          2409817111,
          2500171864,
          2581367074,
          2659444322,
          2722469842,  # chr20
          2781598825,
          2829728720,  # chr22
          2881033286,  # chrX (w/o PAR)
          3033334811,  # chrY (w/o PAR)
          3089739342,  # chrXY (PAR1 & PAR2)
          3095677412)

.hg19.breaks = c(124625310,
                 370850307,
                 591461209,
                 786049562,
                 972084330,
                 1148099493,
                 1313226358,
                 1465977701,
                 1609766427,
                 1748140516,
                 1883411148,
                 2017840353,
                 2142351240,
                 2253610949,
                 2358551415,
                 2454994487,
                 2540769469,
                 2620405698,
                 2690957082,
                 2752034333,
                 2805663772,
                 2855381003,
                 2957184048,
                 3061537076,
                 3092708377,
                 3095685697)

ggchicago <- function (data1, data2, chr = "CHR", bp = "BP", P = "P", logP = TRUE, 
          build = "hg19", significance = c(5e-08), theme_base = theme_publication(), 
          scale_color = scale_color_traditional()) 
{
  requireNamespace("ggplot2")
  requireNamespace("ggrastr")
  require("ggman")
  if (is.null(data1[[bp]]) || is.null(data2[[bp]])) {
    stop("NULL BP")
  }
  if (is.null(data1[[P]]) || is.null(data2[[P]])) {
    stop("NULL P")
  }
  if (is.function(theme_base)) {
    theme_base = theme_base()
  }
  if (is.function(scale_color)) {
    scale_color = scale_color()
  }
  conv1 = .convert2posX(data1[[chr]], data1[[bp]], build)
  data1$x = conv1$posX
  data1$color = as.factor(data1[[chr]])
  data1$y = if (logP) 
    -log10(data1[[P]])
  else data1[[P]]
  conv2 = .convert2posX(data2[[chr]], data2[[bp]], build)
  data2$x = conv2$posX
  data2$color = as.factor(data2[[chr]])
  data2$y = if (logP) 
    log10(data2[[P]])
  else -data1[[P]]
  breaks = if (length(conv1$breaks) > length(conv2$breaks)) 
    conv1$breaks
  else conv2$breaks
  labels = if (length(conv1$labels) > length(conv2$labels)) 
    conv1$labels
  else conv2$labels
  xlabel = ifelse(nchar(conv1$xlabel) > nchar(conv2$xlabel), 
                  conv1$xlabel, conv2$xlabel)
  data = rbind(subset(data1, select = c("x", "y", "color")), 
               subset(data2, select = c("x", "y", "color")))
  plt = ggplot(data, aes(x, y, color = color)) + 
    #ggrastr::rasterise(geom_point(), dpi = 72) + 
    #geom_point()+
    ggrastr::geom_point_rast()+
    geom_hline(yintercept = c(0, -log10(significance),  log10(significance))) + 
    geom_hline(yintercept = c(0, -log10(1e-05), log10(1e-5)), linetype = "dashed") + 
    scale_x_continuous(breaks = breaks, labels = labels) + 
    theme_base + scale_color + theme(legend.position = "none") + 
    xlab(xlabel) + ylab("-log10 P")
  return(plt)
}

#Miami plot
print("Miami")
read_logistic <- function(path){
  df <- fread(cmd = paste('egrep "Z_STAT|ADD"', path),header = T)
  df <- df %>% separate(ID, c('chrpos','SNP'), sep = '%', remove = F)
  df <- df %>% rename('CHR' = '#CHROM', 'BP' = 'POS')
}
LB1 <- read_logistic('/vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme1_500FG.pheno.glm.logistic.hybrid')
LB2 <- read_logistic('/vol/projects/CIIM/Lyme_GWAS/GWAS/assoc/out/lyme2_300BCG.pheno.glm.logistic.hybrid')
pl <- ggchicago(LB2, LB1, scale_color = scale_color_manual(values = rep(c("deepskyblue4", "lightskyblue3"), 11)))+
        theme(axis.text.x = element_text(angle=90))
ggsave("../SVG/Fig1.Miami.svg", plot = pl, width = 14, height = 8)


#qqplot
qqplot <- function(assoc, name){
  chisq <- qchisq(1-na.omit(assoc$P), 1)
  lambda = median(chisq)/qchisq(0.5,1)
  
  print("QQplot")
  png(paste0("../SVG/Fig1.qqplot.",name,".png"))
  qqman::qq(assoc$P, main = paste("lambda=", format(round(lambda, 2), nsmall = 2)))
  dev.off()
}
qqplot(LB1,'lyme1')
qqplot(LB2, 'lyme2')
