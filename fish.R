#!/usr/bin/env Rscript

# ================
# = fish.R =
# ================
#
# loads gene count data and performs differential gene expression
#
# uses a control file, with one analysis per line; columns (tab-separated) are:
# [1] analysis prefix (path to the output folder)
# [2] sample data file (matrix with sample names and information about factors)
# [3] counts file (e.g. STAR_gene_fragment_counts.tsv prepared by STARexplore.R)
# [4] design formula for DESeq2 (defines which factors to use from sample matrix)
# [5] contrast (comparison of interest, comma-separated)
# [6] genes of interest (comma-separated list)

# imports and reading command line arguments
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
control_file <- args[1]

conn <- file(control_file, open = "r")
lines <- readLines(conn)

sample_color_vector <- c('#e31a1c','#ff7f00','#1f78b4') # Vr, VvA, VvI

# main loop, per line of control file
for (i in 1:length(lines)){

    # getting analysis parameters
    params <- strsplit(lines[i], "\t")[[1]]
    print(params)

    prefix <- params[1]
    sample_data_file <- params[2]
    input <- params[3]
    design <- params[4]
    contrast <- params[5]
    genes <- strsplit(params[6], ",")[[1]]

    # reading in data
    sample_data <- read.delim(sample_data_file, row.names= NULL)
    sample_list <- sample_data$sample

    # building count matrix and respective DESeq object
    count_matrix <- read.delim(input)
    count_matrix <- count_matrix[, c("gene", sample_list)] # in case we are not considering all samples in the count matrix
    dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = sample_data, design = as.formula(design), tidy = TRUE)

    # !!! debatable setting #1:
    # pre-filtering rows (genes) with less than 10 reads
    dds <- dds[ rowSums(counts(dds)) >= 10, ]

    # differential expression analysis
    dds <- DESeq(dds)

    # extracting results for defined contrast
    # using alpha=0.05; influences independent filtering (and therefore padj)
    res <- results(dds, contrast = strsplit(contrast, ",")[[1]], alpha = 0.05)

    # writing analysis outputs:

    # create output directory
    dir.create(prefix, recursive = TRUE)

    # export summary(results) as a text file
    sink(file = file.path(prefix, "DESeq2_summary.txt"))
    summary(res)
    sink()

    # export normalized gene counts as csv
    dat <- counts(dds, normalized=TRUE)
    dat <- rownames_to_column(as.data.frame(dat))
    colnames(dat)[1] <- "gene"
    write.csv(as.data.frame(dat), file.path(prefix, "DESeq2_normalized.csv"), quote = FALSE, row.names = FALSE)

    # export differential expression results for all genes as csv
    res2 <- rownames_to_column(as.data.frame(res))
    colnames(res2)[1] <- "gene"
    write.csv(res2, file.path(prefix, "DESeq2_results_all.csv"), quote = FALSE, row.names = FALSE)

    # export differential expression results for genes of interest only
    res3 <- rownames_to_column(as.data.frame(subset(res, rownames(res) %in% genes)))
    colnames(res3)[1] <- "gene"
    write.csv(res3, file.path(prefix, "DESeq2_results_genes.csv"), quote = FALSE, row.names = FALSE)

    # export differential expression results for significant genes only
    res4 <- rownames_to_column(as.data.frame(subset(res, res$padj < 0.05)))
    write.csv(res4, file.path(prefix, "DESeq2_results_padj0.05.csv"), quote = FALSE, row.names = FALSE)

    # print histogram of p-values
    pdf(file.path(prefix, "pvalue_hist.pdf"))
    hist(res$pvalue, breaks = 0:20/20, col = "grey50", border = "white")
    dev.off()

    # print MA plot
    main_factor <- strsplit(contrast, ",")[[1]][1]
    for (coef in resultsNames(dds)[-1]) {
        coef_factor <- strsplit(coef, "_")[[1]][[1]]
        if (main_factor == coef_factor) {
            resLFC <- lfcShrink(dds, coef = coef, type = 'apeglm')
            pdf(file.path(prefix, paste0(coef, "_MA_plot.pdf")))
            plotMA(resLFC, main = coef)
            dev.off()
        }
    }

    # printing count plots for genes of interest
    dir.create(file.path(prefix, "gene_counts_plots"), recursive = TRUE)
    for (gene in genes) {
        tryCatch({

            gene_counts <- plotCounts(dds, gene = gene, intgroup = colnames(colData(dds)), returnData = TRUE)

            # plot based on https://stackoverflow.com/questions/25632242/filled-and-hollow-shapes-where-the-fill-color-the-line-color

            gene_count_plot <- ggplot(gene_counts, aes_string(x = main_factor, y = "count", color = "group", shape = "location")) +
                geom_point(size = 5) +
                scale_y_log10() +
                ggtitle(paste0("normalized counts for gene ", gene, "\n design: ", design)) +
                theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
                scale_color_manual(name = 'Group', labels = c('V. rueppelli','V. vulpes Africa','V. vulpes Iberia'), values = sample_color_vector)

            ggsave(plot = gene_count_plot, file.path(prefix,"gene_counts_plots" , paste0(gene, "_counts.pdf")), width = 22, height = 18, units = "cm", dpi = "retina")

        }, error = function(e){write(paste0("something went wrong for gene ", gene), stdout())})
    }


}
close(conn)
