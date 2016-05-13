# Run the picrust portion of the VMMP pipeline


rule all: 
	input: "picrust/metagenome_predictions_L3.biom"

rule rename_seqs:
	input: "all_seqs.fasta"
	output: "all_seqs_picrust.fasta"
	shell: "cat {input} | sed 's/barcodelabel=//g' | sed 's/;/_/g' > {output}"

rule closed_otu_table:
	input: rules.rename_seqs.output
	output: "picrust/otu_table.biom"
	params: dir="picrust", qiime_params="qiime_params.txt"
	shell: "source activate qiime; \
	        pick_closed_reference_otus.py -i {input} -o {params.dir} -a -O {threads} -p {params.qiime_params}" 

rule normalize:
	input: rules.closed_otu_table.output 
	output: "picrust/normalized_otus.biom"
	shell: "normalize_by_copy_number.py -i {input} -o {output} "

rule predict:
	input: rules.normalize.output
	output: "picrust/metagenome_predictions.biom"
	shell: "predict_metagenomes.py -i {input} -o {output}"

rule categorize:
	input: rules.predict.output
	output: "picrust/metagenome_predictions_L3.biom"
	shell: "categorize_by_function.py -i {input} -c KEGG_Pathways -l 3 -o {output}"
