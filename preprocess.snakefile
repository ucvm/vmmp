
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
