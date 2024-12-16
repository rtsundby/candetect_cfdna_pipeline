#########1#########2#########3#########4#########5#########6#########7#########8
# Mathios, et al. 2021. https://doi.org/10.1038/s41467-021-24994-w
# "In addition to fragmentation profiles, we also computed z-scores for
#  chromosomal arms and a genome-wide summary of the overall cfDNA fragmentation
#  as previously described23,51 with the modifications indicated below. We only
#  analyzed the 39 chromosomal arms that were not acrocentric. Z-scores for each
#  of the 39 autosomal arms were obtained by centering and scaling the total
#  GC-adjusted fragment count for each arm by the mean and standard deviation of
#  the corresponding arm-specific counts in the 54 non-cancer samples used as a
#  reference set."

# For unit testing
args = commandArgs(trailingOnly = TRUE)
cytobands_tsv      = args[1]
libraries_full     = args[2]
healthy_term       = args[3]
counts_tsv         = args[4]
output_csv         = args[5]

# Load required packages and data 
library(tidyverse)
library(data.table)

cytobands = read_tsv(cytobands_tsv, col_names = c("chr","start","end","band","stain"))
counts = read_tsv(counts_tsv)
libraries_full = read_tsv(libraries_full) 

# Get healthy control libraries
healthy_libs =
  libraries_full %>%
  filter(cohort == healthy_term) %>%
  pull(library)
healthy_libs

# Define non-acrocentric arms and their genomic coordinates
(autosomes = c(paste0("chr", seq(1:22))))
(noacro_chroms = c(paste0("chr", seq(1:12)), paste0("chr", seq(16, 20, by =1))))
(acro_chroms = autosomes[!autosomes %in% noacro_chroms])
noacro = cytobands %>%
  mutate(arm = substr(band, 1, 1)) %>%
  group_by(chr, arm) %>%
  summarize(armstart = min(start),
            armend = max(end)) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  pivot_wider(names_from = arm, values_from = c(armstart,armend)) 
noacro

################################################################################
# Perform data transforms specified in lucas workflow:
# https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/394576177f78961816e7adc486d275897be39015/code/preprocessing/04-getZscores.r#L46

# LOESS correction for GC content
# See: https://github.com/cancer-genomics/reproduce_lucas_wflow/blob/394576177f78961816e7adc486d275897be39015/code/PlasmaTools.lucas/R/aneuploidy_zscores_pascores.R#L107
gcCorrectLoess <- function(counts.mult, gc) {
  lower <- 0
  upper <- quantile(counts.mult, .99, na.rm=TRUE)
  trend.points <- counts.mult > lower & counts.mult < upper
  trend.counts <- counts.mult[trend.points]
  trend.gc <- gc[trend.points]
  num.sample.points <- min(length(trend.counts), 10E3L)
  samp <- sample(1:length(trend.counts), num.sample.points)
  #    pad counts and GC
  include <- c(which(gc==min(gc)), which(gc==max(gc)))
  trend.counts <- c(counts.mult[include], trend.counts[samp])
  trend.gc <- c(gc[include], trend.gc[samp])
  initial.trend <- loess(trend.counts ~ trend.gc)
  i <- seq(min(gc, na.rm=TRUE), max(gc, na.rm=TRUE), by = 0.001)
  final.trend <- loess(predict(initial.trend, i) ~ i)
  counts.pred <- predict(final.trend, gc)
  return(counts.pred)
}

# Use non-acrocentric arm positions to count fragments per arm
arm_counts =
  counts %>% 
  # Get total counts, short and long, per position
  group_by(library,chr,start,end,gc) %>% 
  summarize(cov=sum(count)) %>%
  # Annotate arms and filer to non-acro arms
  filter(chr %in% autosomes) %>% left_join(noacro, by = "chr") %>%
  mutate(arm = ifelse(end < armstart_q, "p", "q")) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  subset(select=c(library,chr,start,end,gc,arm,cov))
arm_counts

# Corrections over bins
arm_counts_DT <- setDT(arm_counts)
arm_counts_DT[,counts.mult:=log2((get("cov") + 1))]
arm_counts_DT[,loess.pred:=gcCorrectLoess(counts.mult, gc),by="library"]
arm_counts_DT[,adjusted:=counts.mult - loess.pred]
arm_counts_DT[,adjusted.cent:=adjusted - median(adjusted, na.rm=TRUE), by="library"]
arm_counts_DT[]

# Aggregate bins to arms
arms <- arm_counts_DT[,.(armmean=mean(adjusted.cent, na.rm=TRUE)),
             by=.(library, chr, arm)] 
arms <- arms[,armmean:=armmean - mean(armmean, na.rm = TRUE)]
arms[]

arms_tib <- as_tibble(arms)
arms_tib

# Get mean and sd of healthy libs
healthy_arms =
  arms_tib %>%
  filter(library %in% healthy_libs) %>%
  group_by(chr,arm) %>%
  summarize(mean = mean(armmean, na.rm = T), sd = sd(armmean, na.rm=T)) 
healthy_arms

# z-score against the healthy cohort
arm_z =  
  arms_tib %>%
  left_join(healthy_arms, by = c("chr","arm")) %>%
  mutate(z = (armmean-mean)/sd) %>%
  select(!c("mean", "sd"))
arm_z

arm_z2 =
  arm_z %>%
  select(-armmean) %>%
  unite(chr_arm,chr:arm,sep="_") %>%
  spread(chr_arm,z) %>%
  column_to_rownames('library')
arm_z2

write.csv(arm_z2, file = output_csv, row.names=T)
