#!/usr/bin/env Rscript

# =================
# = STARexplore.R =
# =================
#
# loads relevant STAR results and performs exploratory analyses:
# - stacked bar plot with counts in the different alignment categories
# - table of the same
# - dendrograms and MDS plots based on pairwise Euclidean distances between samples
# - PCA
# - table of gene fragment counts (for importing into Differential Expression scripts)
#
# command line arguments:
# [1] sample data file (tsv with columns: sample, location, wc, raw_reads_M, trimmed_reads_M)
# [2] STAR output folder (containing *_ReadsPerGene.out.tab files)
# [3] output prefix


# 0. imports and reading command line arguments
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
sample_data_file <- args[1]
star_folder <- args[2]
prefix <- args[3]

# 1. setup: reading data and checking
sample_data <- read.delim(sample_data_file, row.names= NULL)
sample_list <- sample_data$sample
files <- paste0(star_folder, sample_list, "_ReadsPerGene.out.tab")

# checking if all files exist, stop script if not
for(file in files){
    if(!file.exists(file)){
        write(paste(file, "doesn't exist!"), stdout())
        stop()
        }
    }

dir.create(prefix)

# reading in gene fragment counts for all samples
# considers the counts in the third column of the STAR file ('forward'/read 1 stranded counts)
temp_list <- list()
for(i in seq_len(length(sample_list))){
    sample <- sample_list[i]
    file <- files[i]
    # colClasses defines which columns to skip ("NULL" in this case skips 2nd and 4th columns)
    temp_d <- read.delim(file, row.names = 1, header = FALSE, colClasses = c(NA, "NULL", NA, "NULL"))
    colnames(temp_d)[1] <- sample
    temp_list[[sample]] <- temp_d
}

counts <- do.call(cbind, temp_list)


# 2. manipulating and exporting fragment counts tables
# separating into different datasets: one for category counts/plots (cat_counts), other for DE (de_counts)
# first 4 rows of the imported data contain categpry counts for N_unmapped, N_multimapping, N_noFeature, N_ambiguous
cat_counts <- counts[seq(1,4),]
de_counts <- counts[-seq(1,4),]
de_counts <- de_counts[!(rownames(de_counts) %in% "transcript_id"),] # cleaning a false transcript called "transcript_id"
# generating N_overlapping by summing all gene overlapping counts
cat_counts <- rbind(cat_counts, N_overlapping = colSums(de_counts))
cat_counts2 <- rownames_to_column(as.data.frame(t(cat_counts)), "sample")

# exporting fragment alignment category counts
write_tsv(cat_counts2, file.path(prefix, "STAR_category_counts.tsv"))

# exporting fragment counts per gene
write_tsv(rownames_to_column(de_counts, "gene"), file.path(prefix, "STAR_gene_fragment_counts.tsv") )


# 3. fragment alignment category plot
cat_counts2$sample <- factor(cat_counts2$sample)
cat_counts2_long <- pivot_longer(cat_counts2, names_to = "category", values_to = "count", N_unmapped:N_overlapping)
cat_counts2_long$count <- cat_counts2_long$count/1000000
cat_counts2_long$category <- as.factor(cat_counts2_long$category)
cat_counts2_long <- cat_counts2_long %>% mutate(category = fct_relevel(category, "N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous", "N_overlapping"))
cat_plot <- ggplot(data=cat_counts2_long, aes(x=sample, y=count, fill=category)) +
    geom_col(position="stack") +
    ylab("number of read pairs (millions)") +
    scale_fill_brewer(name = "STAR count category", labels = c("unmapped","multimapped","no feature","ambiguous feature","overlapping genes"), palette = "RdYlBu")

ggsave(plot = cat_plot, file.path(prefix, "STAR_category_counts.pdf"), dpi = "retina")


# 4. importing & manipulating counts into DESeq2 objects
STAR_se <- SummarizedExperiment(assays = list(STAR_counts = as.matrix(de_counts)), rowData = DataFrame(gene_id=rownames(de_counts)), colData = sample_data)
star_dds <- DESeqDataSet(STAR_se, design = ~ 1)

# transforming counts with rlog
star_rld <- rlog(star_dds)


# 5. Euclidean distance plots (dendrograms and MDS)
# calculating Euclidean distance between pairs of samples
star_dists <- dist(t(assay(star_rld)))

# dendrogram
pdf(file.path(prefix, "STAR_EuclDist_dendrogram.pdf"))
plot(hclust(as.dist(star_dists)), main='', sub='', xlab='', ylab='Euclidean distance')
dev.off()

# MDS plot
star_mds <- cbind(cmdscale(star_dists), as.data.frame(colData(star_dds)))

mds_plot <- ggplot(star_mds, aes(`1`, `2`, color = species, shape = location)) +
    xlab('') +
    ylab('') +
    geom_point(size = 8, aes(fill = species, alpha = factor(ifelse(wc == "captive", 1, 0)))) +
    geom_point(size = 8) +
    scale_shape_manual(values = c(21, 24)) +
    scale_alpha_manual(name = "wc", values = c("1"=0, "0"=1), labels = c("wild", "captive")) +
    geom_label_repel(label = rownames(star_mds))

ggsave(plot = mds_plot, file.path(prefix, "STAR_EuclDist_MDS.pdf"), width = 18, height = 18, units = "cm", dpi = "retina")


# 6. PCA
star_pcadata <- plotPCA(star_rld, intgroup = colnames(colData(star_dds)), returnData = TRUE )
star_percentVar <- round(100 * attr(star_pcadata, "percentVar"))

pca_plot <- ggplot(star_pcadata, aes(x = PC1, y = PC2, color = species, shape = location)) +
    geom_point(size = 8, aes(fill = species, alpha = factor(ifelse(wc == "captive", 1, 0)))) +
    geom_point(size = 8) +
    scale_shape_manual(values = c(21, 24)) +
    scale_alpha_manual(name = "wc", values = c("1"=0, "0"=1), labels = c("wild", "captive")) +
    xlab(paste0("PC1: ", star_percentVar[1],"% variance")) +
    ylab(paste0("PC2: ", star_percentVar[2],"% variance")) +
    geom_label_repel(label = rownames(star_pcadata))

ggsave(plot = pca_plot, file.path(prefix, "STAR_PCA_PC1-PC2.pdf"), width = 18, height = 18, units = "cm", dpi = "retina")
