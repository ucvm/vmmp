# *******************************
# * Vet med microbiome pipeline *
# *******************************

# **** Modules ****

# if these are in the path already then comment this section out
exec(open("/usr/share/Modules/init/python3.py").read())

module('load', 'applications/usearch/9.2.64')
module('load', 'applications/ssu-align/0.1.1')

# **** Variables ****

configfile: "config.yaml"

VERSION='1.3.0'

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


# **** Includes ****

include: "preprocess.snakefile"
include: "{type}.snakefile".format(type = config["analysis_type"])
include: "taxonomy.snakefile"
include: "misc.snakefile"


# **** Rules ****

rule all:
    input: "results/phyloseq.rds", "stats/quality_stat.txt", "results/otu_table.txt",
           "results/QC_report_{}.html".format(config["run_name"])




