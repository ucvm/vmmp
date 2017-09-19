
rule make_otutab:
    input: seqs = rules.merge_pairs.output.merged, db = rules.make_otus.output
    output: biom = "results/otu_table.biom", txt = "results/otu_table.txt"
    threads: config["num_threads"]
    shell: "usearch -usearch_global {input.seqs} -db {input.db} -threads {threads} -strand plus -id 0.97 \
             -biomout {output.biom} -otutabout {output.txt}"

rule ssu_align:
    input: rules.make_otus.output
    output: "ssu_out/ssu_out.bacteria.stk"
    params: dir="ssu_out"
    log: "logs/align.log"
    shell: "ssu-align -f {input} {params.dir} &>> {log}"

rule ssu_mask:
    input: rules.make_otus.output
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
    input: seqs=rules.make_otus.output, 
           otutab=rules.make_otutab.output.biom,
           tree=rules.tree.output
    output: "results/phyloseq.rds"
    params: database=config["database"], 
            revcomp=config["revcomp"],
            taxa_type=config["taxa_type"]
    shell: "Rscript --vanilla assign_taxonomy.R {params.database} {input.seqs} {params.revcomp} {input.otutab} {input.tree} {params.taxa_type} {output}"


