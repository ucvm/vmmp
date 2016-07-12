# *******************************
# * Vet med microbiome pipeline *
# *******************************

# **** Modules ****

# if these are in the path already then comment this section out
exec(open("/usr/share/Modules/init/python3.py").read())

module('load', 'applications/usearch/8.1.1861')
module('load', 'applications/ssu-align/0.1.1')
module('load', 'applications/FastTree/2.1.8')


# **** Variables ****

configfile: "config.yaml"

VERSION='1.0.0'

# **** Imports ****

import glob
from snakemake.utils import R, report

# **** Rules ****

rule all:
    input: "results/phyloseq.rds", "stats/quality_stat.txt" 

rule clip_primers:
    input: r1 = lambda wildcards: glob.glob("{directory}/{sample}_*R1*.fastq*".format(directory=config["read_directory"], sample=wildcards.sample)),
           r2 = lambda wildcards: glob.glob("{directory}/{sample}_*R2*.fastq*".format(directory=config["read_directory"], sample=wildcards.sample))
    output: r1="clipped/{sample}_R1.cut", r2="clipped/{sample}_R2.cut"
    log: "logs/{sample}.log"
    shell: "cutadapt -e 0 -O 17 -g {config[fwd_primer]} -G {config[rev_primer]} -a {config[rev_primer_rc]} -A {config[fwd_primer_rc]} -o {output.r1} -p {output.r2} {input.r1} {input.r2} >> {log}"

rule merge_pairs:
    input: expand("clipped/{sample}_R1.cut", sample=config["samples"]) 
    output: "processed/merged.fq" 
    log: "logs/merged.log"
    threads: config["num_threads"]
    shell: "usearch -fastq_mergepairs {input} -fastqout {output[0]} -relabel @ -threads {threads} -report {log}"

rule calc_stats:
    input: rules.merge_pairs.output
    output: "stats/quality_stat.txt"
    threads: config["num_threads"]
    shell: "usearch -fastq_eestats2 {input} -output {output} -threads {threads} -length_cutoffs '100,*,10'"

rule filter:
    input: rules.merge_pairs.output
    output: "processed/filtered.fa"
    log: "logs/filter.log"
    threads: config["num_threads"]
    shell: "usearch -fastq_filter {input} -fastaout {output} -fastq_maxee {config[max_expected_error]} -fastq_trunclen {config[read_length_cutoff]} -threads {threads} -log {log}"

rule dereplicate:
    input: rules.filter.output
    output: "cluster/uniques.fa"
    threads: config["num_threads"]
    shell: "usearch -derep_fulllength {input} -relabel Uniq -fastaout {output} -sizeout"

rule cluster:
    input: rules.dereplicate.output
    output: "results/otus.fa"
    log: "logs/cluster.log"
    shell: "usearch -cluster_otus {input} -otus {output} -relabel OTU -minsize 2 -log {log}"

rule make_otutab:
    input: seqs=rules.merge_pairs.output, db=rules.cluster.output
    output: "results/otu_table.biom"
    threads: config["num_threads"]
    shell: "usearch -usearch_global {input.seqs} -db {input.db} -threads {threads} -strand plus -id 0.97 -biomout {output}"

rule ssu_align:
    input: rules.cluster.output
    output: "ssu_out/ssu_out.bacteria.stk"
    params: dir="ssu_out"
    log: "logs/align.log"
    shell: "ssu-align -f {input} {params.dir} &>> {log}"

rule ssu_mask:
    input: rules.ssu_align.output
    output: "ssu_out/ssu_out.bacteria.mask.afa"
    params: dir="ssu_out"
    log: "logs/align.log"
    shell: "ssu-mask --afa {params.dir} &>> {log}"

rule tree:
    input: rules.ssu_mask.output
    output: "results/otus.tre"
    log: "logs/tree.log"
    shell: "FastTree -nt {input} >{output} &> {log}"

rule assign_taxonomy:
    input: seqs=rules.cluster.output, 
           otutab=rules.make_otutab.output,
           tree=rules.tree.output
    output: "results/phyloseq.rds"
    params: database=config["database"], 
            revcomp=config["revcomp"],
            taxa_type=config["taxa_type"]
    run:
        R('''
            library(dada2)
            library(phangorn)
            library(ape)
            library(phyloseq)
            library(stringr)
            library(Biostrings)

            database = "{params.database}"
            otu_seqs = "{input.seqs}"
            revcomp = {params.revcomp}
            otu_tab = "{input.otutab}"
            tree_file = "{input.tree}"
            taxa_type = "{params.taxa_type}"
            out_file = "{output}"

            tree = read.tree(tree_file)
            tree = midpoint(tree)
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

            tax_table(data) = tax_table(taxa)

            saveRDS(data, out_file)
            
        ''')

rule clean:
    shell: "rm -rf clipped cluster logs processed results ssu_out stats"

rule print_pipeline_code:
    run: 
        shell("echo '# --------------------------\n# Vet Med Microbiome Pipeline\n# -----------------------\n'")
        shell("echo '# Version: {VERSION}'")
        shell("echo -n '# Date: '; date +'%m-%d-%Y'")
        shell("echo '\n# Versions:\n# ---------'")
        shell("echo -n '# '; usearch --version")
        shell("echo -n '# cutadapt '; cutadapt --version")
        shell("echo -n '# '; FastTree 2>&1 | grep version | sed -r 's/Usage for //' ")
        shell("ssu-align -h 2>&1 | grep SSU-ALIGN")
        shell("echo '\n# Code:\n'")
        shell("snakemake -n -p -q")


