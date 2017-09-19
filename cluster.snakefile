
rule make_otus:
    input: rules.dereplicate.output
    output: "results/otus.fa"
    log: "logs/cluster.log"
    threads: config["num_threads"]
    shell: "usearch -cluster_otus {input} -otus {output} -relabel OTU -minsize 2 -log {log} -threads {threads}"
