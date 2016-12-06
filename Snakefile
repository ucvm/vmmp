# *******************************
# * Vet med microbiome pipeline *
# *******************************

# **** Modules ****

# if these are in the path already then comment this section out
exec(open("/usr/share/Modules/init/python3.py").read())

module('load', 'applications/usearch/9.0.2132')
module('load', 'applications/ssu-align/0.1.1')
#module('load', 'applications/FastTree/2.1.8')
#module('load', 'applications/fastqc/0.11.3')

# **** Variables ****

configfile: "config.yaml"

VERSION='1.2.0'

# **** Conda ****

# set location of conda install - being explicit about this avoids problems with multiple conda installs
CONDA="~/miniconda3/bin/"

# assumes env called 'vmmp' already created, use the provided vmmp_env.yml file to create
ENV="source {CONDA}/activate vmmp".format(CONDA = CONDA)

# prefix shell commands with the above source command
shell.prefix(ENV + '; ')


# **** Imports ****


import glob
from snakemake.utils import R, report


# **** Rules ****

rule all:
    input: "results/phyloseq.rds", "stats/quality_stat.txt", "results/otu_table.txt",
           "results/QC_report_{}.html".format(config["run_name"])

rule clip_primers:
    input: r1 = lambda wc: glob.glob("{dir}/{sample}_*R1*.fastq*".format(dir=config["read_directory"], sample=wc.sample)),
           r2 = lambda wc: glob.glob("{dir}/{sample}_*R2*.fastq*".format(dir=config["read_directory"], sample=wc.sample))
    output: r1="clipped/{sample}_R1.cut", r2="clipped/{sample}_R2.cut"
    log: "logs/{sample}.log"
    shell: """ 
	    cutadapt -e 0 -O 17 -g {config[fwd_primer]} -G {config[rev_primer]} \
            -a {config[rev_primer_rc]} -A {config[fwd_primer_rc]} \
            -m 50 -q {config[q_trim]} \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} >> {log}
            """
rule fastqc:
    input: lambda wc: glob.glob("{dir}/*.fastq*".format(dir = config["read_directory"]))
    output: touch("fastqc.done")
    threads: config["num_threads"] 
    shell: "mkdir -p fastqc; fastqc -t {threads} -o fastqc {input}"

rule multiqc:
    input: "fastqc.done"
    output: "multiqc_data/multiqc_fastqc.txt"
    shell: "multiqc -f fastqc"

rule merge_pairs:
    input: expand("clipped/{sample}_R1.cut", sample=config["samples"]) 
    output: merged = "processed/merged.fq", log = "logs/merged.log"
    threads: config["num_threads"]
    shell: "usearch -fastq_mergepairs {input} -fastqout {output.merged} -relabel @ -threads {threads} -report {output.log}"

rule calc_stats:
    input: rules.merge_pairs.output.merged
    output: "stats/quality_stat.txt"
    threads: config["num_threads"]
    shell: "usearch -fastq_eestats2 {input} -output {output} -threads {threads} -length_cutoffs '100,*,10'"

rule filter:
    input: rules.merge_pairs.output.merged
    output: filtered = "processed/filtered.fa", log = "logs/filter.log"
    threads: config["num_threads"]
    shell: "usearch -fastq_filter {input} -fastaout {output.filtered} -fastq_maxee {config[max_expected_error]} -fastq_trunclen {config[read_length_cutoff]} -threads {threads} -log {output.log}"

rule dereplicate:
    input: rules.filter.output.filtered
    output: "cluster/uniques.fa"
    threads: config["num_threads"]
    shell: "usearch -derep_fulllength {input} -relabel Uniq -fastaout {output} -sizeout -threads {threads}"

rule cluster:
    input: rules.dereplicate.output
    output: "results/otus.fa"
    log: "logs/cluster.log"
    threads: config["num_threads"]
    shell: "usearch -cluster_otus {input} -otus {output} -relabel OTU -minsize 2 -log {log} -threads {threads}"

rule make_otutab:
    input: seqs = rules.merge_pairs.output.merged, db = rules.cluster.output
    output: biom = "results/otu_table.biom", txt = "results/otu_table.txt"
    threads: config["num_threads"]
    shell: "usearch -usearch_global {input.seqs} -db {input.db} -threads {threads} -strand plus -id 0.97 \
             -biomout {output.biom} -otutabout {output.txt}"

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
    shell: "FastTree -nt {input} >{output} 2> {log}"

rule assign_taxonomy:
    input: seqs=rules.cluster.output, 
           otutab=rules.make_otutab.output.biom,
           tree=rules.tree.output
    output: "results/phyloseq.rds"
    params: database=config["database"], 
            revcomp=config["revcomp"],
            taxa_type=config["taxa_type"]
    shell: "Rscript --vanilla assign_taxonomy.R {params.database} {input.seqs} {params.revcomp} {input.otutab} {input.tree} {params.taxa_type} {output}"

rule qc_report:
    input: fastqc_data = rules.multiqc.output,
       	   merge_log = rules.merge_pairs.output.log,
	   filter_log = rules.filter.output.log
    output: "results/QC_report_{}.html".format(config["run_name"])
    params: runID = config["run_name"]
    shell: "Rscript --vanilla compile_QC_report.R {output} {input.fastqc_data} {params.runID} {input.merge_log} {input.filter_log}"

rule clean:
    shell: "rm -rf clipped cluster logs processed results ssu_out stats fastqc.done"

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
	shell("conda list")
        shell("echo '\n# Code:\n'")
        shell("snakemake -n -p -q")

rule print_conda_env:
    run:
        shell("conda env export")

