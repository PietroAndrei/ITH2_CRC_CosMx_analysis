# Analysis of CRC CosMx 1K samples
The analysis and data contained in this repository have been originally produced on a Linux HPC within conda environments, using conda v22.9.0 (updated from v4.12.0) and mamba v0.26.0 to manage package dependencies and Snakemake v7.14.0 to run the code and manage the analysis workflow
We recommend to download the same conda version from https://repo.anaconda.com/miniconda/, then follow the instructions on https://docs.anaconda.com/miniconda/ ('Quick command line install' section) according to your operating system

NB** The analyses have been always run in a Linux environment. Therefore, we don't know whether package dependencies might create issues on MacOS or Windows

## 1) Reproduce original environments
Run the following command to recreate the same conda environment used for the analysis: 

```console
mamba env create -f code/envs/cosmx_archive.yaml
```

You will need to install Seurat v5.0.0 and InSituType v1.0.0 separately with R:

```r
remotes::install_version("Seurat", version = "5.0.0")
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType@v1.0.0")
```

Create a separate conda environment that will be used to run snakemake (v7.14.0):

```console
mamba env create -f code/envs/snakemake.yaml
```

## 2) Reproduce paper results
The whole pipeline is expected to be run within the `code/` subdirectory:

```console
cd code/ #Transfer within the code/ subdirectory from the Figure4/ directory
ln -s snakefiles/Snakefile Snakefile #Create a soft symlink to the snakefiles/Snakefile file in the current directory
conda activate snakemake #Activate conda environment with Snakemake installed in it
snakemake -np # Execute a Snakemake (-np) dry run to check for any bug/error present within the Snakefile
snakemake --use-conda # Execute Snakemake by telling it to activate the conda environment specified within each rule (ie 'cosmx') before running any script
```

***Some steps of the Snakemake pipeline require a considerable amount of memory. Therefore, it is recommended to run the pipeline on a HPC environment. Alternatively, you can run specific pipeline rules by explicitly call their name***

```console
snakemake -np rule_name
snakemake rule_name
```

## Practical suggestions for HPC users
If you have access to a HPC environment to run the pipeline (recommended), Snakemake can directly communicate with job schedulers like SLURM and SGE using Snakemake profiles (https://github.com/Snakemake-Profiles). 
You can download our modified version of the Snakemake profile for SGE scheduler at https://github.com/PietroAndrei/SnakemakeProfile_CSSGE. 
Follow the instructions on the Github repository to download and set up the Snakemake profile on the cluster. 
Once the Snakemake profile is accessible in the right `~/.config/snakemake/` folder, you can run Snakemake using the following command:

```console
snakemake --profile cssge -j 1 --use-conda
```

In the previous code, `--profile` is used to specify the name of the Snakemake profile that will be used to manage snakemake jobs (the name must match with the profile folder name in the `~/.config/snakemake/` directory).
`-j` or `--jobs` argument specifies the number jobs that will be submitted in parallel on the cluster. 
Memory and runtime requirements for job submission are specified within each rule of the Snakefile in the "resources" section.
