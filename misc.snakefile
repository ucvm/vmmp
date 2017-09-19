
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
