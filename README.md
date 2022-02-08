# assembly_pipeline

The script assembly_pipeline.py is a computational pipeline to perform _de novo_ assembly of bacterial genomes from Illumina paired-end reads. The pipeline is based on the assember [SPAdes](https://github.com/ablab/spades) and the [improve_assembly](https://github.com/sanger-pathogens/assembly_improvement) pipeline designed to polish the SPAdes assembly by scaffolding and gap filling.

# Docker Installation

The easiest and recommended way to install and run assembly_pipeline.py is via its Docker implementation.

The Docker image is available on: https://hub.docker.com/r/francesccoll/assembly_pipeline/

# Local Installation

assembly_pipeline.py is Python script that would work provided that all required dependencies below (both python modules and software) are installed in your local machine.

## Required dependencies

### Software
* [fastqcheck](https://github.com/VertebrateResequencing/fastqcheck) version >= 1.1
* [spades.py](https://github.com/ablab/spades) version >= v3.15.3
* [improve_assembly](https://github.com/sanger-pathogens/assembly_improvement)
* [quast.py](https://github.com/ablab/quast) version >= v5.1.0rc1


### Python Modules
* [subprocess](https://docs.python.org/3/library/subprocess.html)


# Usage

```console
usage: assembly_pipeline.py [-h] -1 FASTQ1_FILE -2 FASTQ2_FILE -i SAMPLE_ID -r
                            RESULTS_DIR [-d DELETE_TMP] [--version]
                            [-t THREADS] [-s SPADES_DIR] [-m IMPROVED_DIR]

Pipeline for bacterial de novo assembly using Spades and improve_assembly from
paired Illumina data

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -1 FASTQ1_FILE, --forward_reads FASTQ1_FILE
                        fastq file with forward reads
  -2 FASTQ2_FILE, --reverse_reads FASTQ2_FILE
                        fastq file with reverse reads
  -i SAMPLE_ID, --sample_id SAMPLE_ID
                        sample id used as prefix to name output files
  -r RESULTS_DIR, --results_dir RESULTS_DIR
                        directory to store pipeline's final assembly

optional arguments:
  -d DELETE_TMP, --delete_tmp DELETE_TMP
                        delete assembly files (except for contigs.fa)
  --version             show program's version number and exit

spades arguments (optional):
  -t THREADS, --spades_threads THREADS
                        number of threads used by Spades
  -s SPADES_DIR, --spades_dir SPADES_DIR
                        directory to store Spades resulting files
  -m IMPROVED_DIR, --improved_dir IMPROVED_DIR
                        directory to store improve_assembly resulting files
```

# License

assembly_pipeline.py is a free software, licensed under [GNU General Public License v3.0](https://github.com/francesccoll/assembly_pipeline/blob/main/LICENSE)

# Feedback/Issues

Use the [issues page](https://github.com/francesccoll/assembly_pipeline/issues) to report on installation and usage issues.

# Citations
_Not available yet_

