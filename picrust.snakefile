# Run the picrust portion of the VMMP pipeline


rule all: 
	input: "picrust/metagenome_predictions_ko_L3.biom",
               expand("picrust/metagenome_predictions_{type}.biom", type = ["ko", "cog", "rfam"])

rule closed_otu_table:
	input: "processed/filtered.fa"
	output: "picrust/otu_table.biom"
	params: dir="picrust", qiime_params="qiime_params.txt"
	shell: "source activate qiime; \
	        pick_closed_reference_otus.py -i {input} -o {params.dir} -a -O {threads} -p {params.qiime_params} -f" 

rule normalize:
	input: rules.closed_otu_table.output 
	output: "picrust/normalized_otus.biom"
	shell: "source activate picrust; normalize_by_copy_number.py -i {input} -o {output} "

rule predict_ko:
	input: rules.normalize.output
	output: "picrust/metagenome_predictions_ko.biom"
	shell: "source activate picrust; \
                predict_metagenomes.py -i {input} -o {output} -t ko -a ko_nsti_scores.txt "

rule predict_cog:
	input: rules.normalize.output
	output: "picrust/metagenome_predictions_cog.biom"
	shell: "source activate picrust; \
                predict_metagenomes.py -i {input} -o {output} -t cog -a ko_nsti_scores.txt --with_confidence"

rule predict_rfam:
	input: rules.normalize.output
	output: "picrust/metagenome_predictions_rfam.biom"
	shell: "source activate picrust; \
                predict_metagenomes.py -i {input} -o {output} -t rfam -a ko_nsti_scores.txt --with_confidence"

rule categorize_ko:
	input: rules.predict_ko.output
	output: "picrust/metagenome_predictions_ko_L3.biom"
	shell: "source activate picrust; categorize_by_function.py -i {input} -c KEGG_Pathways -l 3 -o {output}"


rule clean:
	shell: "rm -rf picrust/"
