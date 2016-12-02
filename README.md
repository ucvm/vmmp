[![Snakemake](https://img.shields.io/badge/snakemake-≥3.4.1-brightgreen.svg?style=flat-square)](https://bitbucket.org/johanneskoester/snakemake)

# Vet Med Microbiome Pipeline

This is a Snakemake pipeline designed to process 16S gene survey data using the [UPARSE](http://drive5.com/uparse/) OTU clustering method.  It assigns taxonomy to the representative OTUs with the RDP classifer, aligns the sequences with ssu-align and constructs a tree with FastTree.  The pipeline has been developed for the Faculty of [Veterinary Medicine](http://www.vet.ucalgary.ca) at the University of Calgary.

This repository is provided for reference purposes for publications that use this pipeline and is not provided as a tool for others to use.  This means there is no support or help provided.  That being said anyone is welcome to clone the repository and use the pipeline or feel free to use it as a guide to write your own Snakemake pipeline.

## Install

Clone this repository to a location of your choosing.

```
git clone https://github.com/ucvm/vmmp
```

The recommended way to run the pipeline is to use the [conda](https://anaconda.org)  See the snakemake [webpage](https://bitbucket.org/johanneskoester/snakemake/wiki/Home) for details on how to do this with your snakemake install.

Almost all dependencies can be installed using the provided conda environment file `vmmp_env.yml`.  It will create a new environment called `vmmp`.

```
conda env create -f vmmp_env.yml
```

As of now the R package `dada2` needs to be installed directly from Biocondutor due to dependy issues.

Once everything is installed you can: `source activate vmmp` then run snakemake.  Only usearch and ssu-align will need to be installed manually.

## Dependencies

As configured the snakefile will load the required dependencies using environment modules installed on our local server.  As long as the dependencies below are in your path then there is no need to use the modules.  Simply comment those lines out.  Also, you'll need to comment out the `onsuccess` and `onerror` portions or replace with your own code.  The `push` command is custom script to push a notification to my [Pushover](https://pushover.net) account.


- Python: 3 and above
- [Snakemake](https://bitbucket.org/johanneskoester/snakemake): 3.9.0
- [usearch](http://www.drive5.com/usearch/download.html): 8.1.1861
- [cutadapt](http://cutadapt.readthedocs.org/en/stable/): 1.8.3
- R: 3.3.2 with the following packages installed: phangorn, ape, phyloseq, dada2, stringr, Biostrings
- [ssu-align](http://eddylab.org/software/ssu-align/): 0.1.1
- [FastTree](http://www.microbesonline.org/fasttree/): 2.1.8

### A note on the taxonomy databases

You'll need a copy of your database of choice formatted to be used by dada2::assignTaxonomy.  You can make this yourself or use one provided by the dada2 authors (see the [documentation](http://benjjneb.github.io/dada2/assign.html)).

## Config file

The pipeline requires a config file, written in yaml, to run.  See the provided example file. Most options are self-explanatory and simple to setup.  Example primer sets are given for common protocols at our institution - these can be changed as required.  Sample names should be unique and contained within the file name.

## Quality check

As of now the pipeline requires manual inspection of the quality data to determine the best parameters for quality filtering.  This is done by filtering a single sample with a range of different parameters and inspecting the results to determine the optimal setting for the expected error (-fastq_maxee) and truncation length (-fastq_trunclen) parameters provided to the usearch -fastq_filter command.  

The quality check is run with `snakemake calc_stats` which runs the pipeline up to `calc_stats` rule. The `quality_stats.txt` file in the stats folder will contain the results.

## Running

Once the quality filtering parameters have been determined and the config file constructed the pipeline can be tested with `snakemake -n -p` which will print out the commands to be run without actually running them.  If all looks good you can run the pipeline with `snakmake` or add the `-j` option with the required number of cores.  If you want to run the pipeline on your local cluster you can do that too as snakemake has cluster support built in (see the snakemake documentation). 


## Results

There are various intermediate folders including a folder with log files that can be inspected if an error is encountered.  The main output is in the 'results' folder. The phyloseq.rds file is an R loadable file that contains a phyloseq object ready to analyze with the otu table, OTU sequences, taxonomy, and phylogenetic tree all pre-loaded.


## Pipeline summary

Most of the preprocessing steps for creating the OTU table are as outlined on the UPARSE [webpage](http://drive5.com/usearch/manual/uparse_pipeline.html).  The basic steps are as follows.

1. Clipping the forward and reverse 16S primers, and any adaptor contamination, with cutadapt
2. Merge the forward and reverse reads with usearch
3. Filter with expected error method and truncate sequences at fixed length
5. Dereplicate with usearch
6. Cluster with `usearch -cluster_otus -minsize 2`
8. Map reads to OTUs with `usearch -usearch_global -biomout`
9. Align OTUs with ssu-align and mask with ssu-mask
10. Build tree with FastTree 
11. Assign taxonomy with RDP classifer as implemented in `dada2::assignTaxonomy`, using the specified database 
12. Load all results into phyloseq object ready for analysis 

## Provenance

To get a list of all the versions of the software used along with the pipeline version and a list of shell commands run by the pipeline type `snakemake print_pipeline_code`.

## Future development

This pipeline will evolve as the analysis tools for 16S data evolve.  New tools and features will be developed in a separate branch, with master remaining stable.  

### Picrust

Support is being added to generate a [PICRUSt](http://picrust.github.io/picrust/) analysis.  This is picrust.Snakefile and it takes the filtered and merged reads from the main pipeline to create a 'closed reference' OTU table with Greengenes as the reference.  This is the only way to run picrust (as per their documentation) and although potentially useful will need to be interpreted carefully.

Picrust analysis depends on Qiime 1.9.1 and PICRUSt 1.0.0

