
rule make_otus:
    input: rules.dereplicate.output
    output: "results/otus.fa"
    log: "logs/cluster.log"
    threads: config["num_threads"]
    shell: "usearch -unoise2 {input} -fastaout {output} -threads {threads}"
    
