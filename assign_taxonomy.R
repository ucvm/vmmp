#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(dada2)
library(ape)
library(phyloseq)
library(stringr)
library(Biostrings)
library(readr)
library(dplyr)

# install phangorn if needed because conda doesn't have it yet
if (!require("phangorn")) install.packages("phangorn", repos = "https://cran.rstudio.com")

database = args[1]
otu_seqs = args[2]
revcomp = args[3] 
otu_tab = args[4] 
tree_file = args[5] 
taxa_type = args[6]
out_file = args[7]

tree = read.tree(tree_file)
tree = phangorn::midpoint(tree)
write.tree(tree, tree_file)

data = import_biom(otu_tab, tree_file, otu_seqs)

seqs = refseq(data)
if (revcomp) seqs = reverseComplement(seqs)
seqs = as.character(seqs)

if (taxa_type == "rdp" | taxa_type == "silva") {{
    taxa = assignTaxonomy(seqs, database)
    rownames(taxa) = names(seqs)
    colnames(taxa) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
}} else if (taxa_type == "gg") {{
    taxa = assignTaxonomy(seqs, database)
    taxa = apply(taxa, 2, function(x) {{
        y = str_match(x, "[kpcofgs]__(.*)")[,2]
        replace(y, y == "", NA)
        }})
    colnames(taxa) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    rownames(taxa) = names(seqs)
}} else {{
    stop("Could not detect taxonomy type.  Should be one of c('rdp', 'silva', 'gg')")
}}

taxa %>% as.data.frame() %>% tibble::rownames_to_column("OTU") %>% 
	write_tsv("results/tax_table.txt")

tax_table(data) = tax_table(taxa)

saveRDS(data, out_file)
