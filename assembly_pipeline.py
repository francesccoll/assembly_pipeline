#!/usr/bin/env python3

import string, re
import argparse
import logging
import os
import subprocess
import sys

# ------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------

_DEPENDENCIES = ['fastqcheck', 'spades.py', 'improve_assembly', 'quast.py']

# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------


def parse_arguments():
    description = "Pipeline for bacterial de novo assembly using Spades and improve_assembly from paired Illumina data"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-1", "--forward_reads", action="store", dest="forward_reads",
        help="fastq file with forward reads", required=True, metavar="FASTQ1_FILE"
    )
    group.add_argument(
        "-2", "--reverse_reads", action="store", dest="reverse_reads",
        help="fastq file with reverse reads", required=True, metavar="FASTQ2_FILE")
    group.add_argument(
        "-i", "--sample_id", action="store", dest="sample_id",
        help="sample id used as prefix to name output files", required=True, metavar="SAMPLE_ID")
    group.add_argument(
        "-r", "--results_dir", action="store", dest="results_dir",
        help="directory to store pipeline's final assembly", required=True, metavar="RESULTS_DIR")

    group = parser.add_argument_group('optional arguments')
    group.add_argument(
        "-d", "--delete_tmp", action="store", dest="delete_tmp",
        help="delete assembly files (except for contigs.fa)", required=False, default=True, metavar="DELETE_TMP")
    group.add_argument(
        '--version', action='version', version='%(prog)s 1.0')

    group = parser.add_argument_group('spades arguments (optional)')
    group.add_argument(
        "-t", "--spades_threads", action="store", dest="spades_threads",
        help="number of threads used by Spades", required=False, default=16, metavar="THREADS")
    group.add_argument(
        "-s", "--spades_dir", action="store", dest="spades_dir",
        help="directory to store Spades resulting files", required=False, metavar="SPADES_DIR")
    group.add_argument(
        "-m", "--improved_dir", action="store", dest="improved_dir",
        help="directory to store improve_assembly resulting files", required=False, metavar="IMPROVED_DIR")

    return parser.parse_args()


def check_fastq_validity(fastq_file):
    """ This function makes sure a fastq file exists """
    is_valid = False
    if os.path.isfile(fastq_file):
        logging.info(f'{fastq_file} found')
        is_valid = True
    else:
        logging.error(f'Cannot find file {fastq_file}!')
    return is_valid


def get_sample_id_from_fastq_file(fastq1, fastq2):
    """ This function extracts the sample Id from the fastq file names """
    sampleId1 = fastq1.split('/')[-1].replace('_1.fastq.gz', '').strip()
    sampleId2 = fastq2.split('/')[-1].replace('_2.fastq.gz', '').strip()
    if not sampleId1 == sampleId2:
        logging.error(f'Fastq files have a different prefix!')
    return sampleId1


def check_arguments(args):
    """ This function sets arguments not defined by the user """
    if not check_fastq_validity(args.forward_reads):
        sys.exit(-1)
    if not check_fastq_validity(args.reverse_reads):
        sys.exit(-1)

    # If sample Id not defined by user, get from fastq file names
    if args.sample_id is None:
        args.sample_id = get_sample_id_from_fastq_file(args.forward_reads, args.reverse_reads)

    # If final improved assembly found, exit with an error
    improved_assembly_final = os.getcwd() + "/" + args.sample_id + ".spades.improved.fasta"
    if os.path.isfile(improved_assembly_final):
        logging.error(f'Final pipeline file found {improved_assembly_final}! Existing...')
        sys.exit(-1)

    # If spades directory not defined by user, create one using sample_id
    if args.spades_dir is None:
        args.spades_dir = os.getcwd() + "/" + args.sample_id + "_spades"

    # If improve assembly directory not defined by user, create one using sample_id
    if args.improved_dir is None:
        args.improved_dir = os.getcwd() + "/" + args.sample_id + "_improved"

    return args


def check_dependency(executable_name):
    """ Returns true if executable exists, else false """
    found = False
    output = subprocess.check_output(['which', executable_name]).strip()
    if output:
        found = True
    return found


def run_fastqcheck(fastq_file):
    """ This function runs fastqcheck and returns the maximum read length """
    uncompress = subprocess.Popen(
        ['gunzip', '-c', fastq_file],
        stdout=subprocess.PIPE,
    )

    fastqcheck = subprocess.Popen(
        ['fastqcheck'],
        stdin=uncompress.stdout,
        stdout=subprocess.PIPE,
    )

    end_of_pipe = fastqcheck.stdout

    read_length = "NA"

    for line in end_of_pipe:
        if "average" in line.decode('utf-8').strip():
            items = line.decode('utf-8').strip().split(',')

    for item in items:
        if "max" in item:
            read_length = item.replace('max', '').strip()

    return read_length


def create_spades_command(read_length, fastq1, fastq2, spades_out_dir, spades_threads):
    """ This function chooses the right k-mer range for Spades from the read length and outputs spades command line"""
    kmer_range = "auto"
    min_read_length = 30
    max_read_length = 301
    read_length = int(read_length)

    if not min_read_length < read_length <= max_read_length:
        logging.error(f'Read length not supported {read_length}!')
        sys.exit(-1)

    if min_read_length < read_length <= 39:
        kmer_range = "15,21,27"
    elif 40 < read_length <= 69:
        kmer_range = "15,21,33"
    elif 70 < read_length <= 124:
        kmer_range = "21,33,55"
    elif 125 < read_length <= 249:
        kmer_range = "21,33,55,77"
    elif 250 < read_length <= max_read_length:
        kmer_range = "21,33,55,77,99,127"

    spades_command = ["spades.py",
                      "-k", kmer_range,
                      "--careful",
                      "-1", fastq1, "-2", fastq2,
                      "-o", spades_out_dir, "-t", str(spades_threads)]

    return spades_command


def run_spades_command(spades_command):
    subprocess.run(
        spades_command,
    )


def create_improve_assembly_command(fastq1, fastq2, spades_assembly, improve_out_dir):
    if not os.path.isfile(spades_assembly):
        logging.error(f'Spades assembly {spades_assembly} not found!')
        sys.exit(-1)

    improve_assembly_command = ["improve_assembly",
                                "-a", spades_assembly,
                                "-f", fastq1, "-r", fastq2,
                                "-o", improve_out_dir]

    return improve_assembly_command


def run_improve_assembly_command(improve_assembly_command):
    subprocess.run(
        improve_assembly_command,
    )


def run_quast_on_assembly(assembly, output_dir):
    if not os.path.isfile(assembly):
        logging.error(f'Input assembly file {assembly} not found!')
        sys.exit(-1)
        
    subprocess.run(
        ["quast.py", assembly,
         "--fast",
         "--threads", "1",
         "-o", output_dir],
    )


def copy_files(file1, file2):
    subprocess.run(
        ["cp", file1, file2],
    )


def delete_directory(directory):
    subprocess.run(
        ["rm", "-r", directory],
    )


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------


def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Check validity of argments. Set if not defined by the user
    args = check_arguments(args)

    logging.info('Making sure dependencies exist...')
    for dependency in _DEPENDENCIES:
        if check_dependency(dependency):
            logging.info(f'{dependency} is installed!')
        else:
            logging.error(f'{dependency} is NOT installed!')
            sys.exit(-1)

    logging.info(f"Running fastqcheck on {args.forward_reads}")
    read_length1 = run_fastqcheck(args.forward_reads)
    if not read_length1.isnumeric():
        logging.error(f'Read length could not be extracted from {args.forward_reads}!')
        sys.exit(-1)
    else:
        logging.info(f"The maximum read length from {args.forward_reads} is {read_length1}")

    logging.info(f"Running fastqcheck on {args.reverse_reads}")
    read_length2 = run_fastqcheck(args.reverse_reads)
    if not read_length2.isnumeric():
        logging.error(f'Read length length could not be extracted from {args.reverse_reads}!')
        sys.exit(-1)
    else:
        logging.info(f"The maximum read length from {args.reverse_reads} is {read_length2}")

    spades_assembly = args.spades_dir + "/" + "contigs.fasta"

    if not os.path.isfile(spades_assembly):
        spades_command = create_spades_command(
            read_length1,
            args.forward_reads,
            args.reverse_reads,
            args.spades_dir,
            args.spades_threads)
        logging.info(f"Running Spades: {spades_command}")
        run_spades_command(spades_command)
    else:
        logging.info(f"Spades assembly {spades_assembly} already found")

    quast_assembly_output = args.spades_dir + '/' + 'transposed_report.tsv'
    if not os.path.isfile(quast_assembly_output):
        logging.info(f"Running quast on : {spades_assembly}")
        run_quast_on_assembly(spades_assembly, args.spades_dir)
    else:
        logging.info(f"Quast output {quast_assembly_output} already found")

    if not os.path.exists(args.improved_dir):
        os.mkdir(args.improved_dir)

    improved_assembly = args.improved_dir + "/" + "scaffolds.scaffolded.gapfilled.length_filtered.sorted.fa"

    if not os.path.isfile(improved_assembly):
        improve_assembly_command = create_improve_assembly_command(
            args.forward_reads,
            args.reverse_reads,
            spades_assembly,
            args.improved_dir)
        logging.info(f"Running improve_assembly: {improve_assembly_command}")
        run_improve_assembly_command(improve_assembly_command)
    else:
        logging.info(f"Improved assembly {spades_assembly} already found")

    quast_improved_output = args.improved_dir + '/' + 'transposed_report.tsv'
    if not os.path.isfile(quast_improved_output):
        logging.info(f"Running quast on : {improved_assembly}")
        run_quast_on_assembly(improved_assembly, args.improved_dir)
    else:
        logging.info(f"Quast output {quast_improved_output} already found")

    # Coping final files
    spades_assembly_final = args.results_dir + os.getcwd() + "/" + args.sample_id + ".spades.contigs.fasta"
    improved_assembly_final = args.results_dir + "/" + args.sample_id + ".spades.improved.fasta"
    quast_assembly_final = args.results_dir + "/" + args.sample_id + ".spades.contigs.quast.csv"
    quast_improved_final = args.results_dir + "/" + args.sample_id + ".spades.improved.quast.csv"

    logging.info(f'Creating final assembly files {spades_assembly_final} and {improved_assembly_final}')
    copy_files(spades_assembly, spades_assembly_final)
    copy_files(improved_assembly, improved_assembly_final)
    copy_files(quast_assembly_output, quast_assembly_final)
    copy_files(quast_improved_output, quast_improved_final)

    # Deleting temporary assembly directories
    if args.delete_tmp:
        logging.info(f'Deleting temporary assembly directories {args.spades_dir} and {args.improved_dir}')
        delete_directory(args.spades_dir)
        delete_directory(args.improved_dir)


if __name__ == "__main__":
    _main()
