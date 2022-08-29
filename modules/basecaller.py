"""
Basecaller wrapper module used to create basecaller commands for bonito and chiron
"""


#Basecaller wrappers for Bonito and Chiron
from Bio import SeqIO

import os
import subprocess
from natsort import natsorted


def bonito(directory, result_directory, outfile='output.fastq', env='bonito', batchsize=0, chunksize=0,
           model='dna_r9.4.1_e8.1_sup@v3.3'):
    """
    --DESCRIPTION--
    Creates a command line that activates the relevant conda environment and runs Bonito
    Conda environment with installed Bonito packages is required before running
    Bonito will produce a single FASTQ file as output

    --PARAMETERS--
    directory: directory path where all files to be basecalled are stored
    result_directory: result directory path - final basecalled sequence is stored here
    outfile: name of the basecalled file
    env: name of the conda environment used to run Bonito
    batchsize: used to control amount of GPU memory used
    chunksize: used to control accuracy of shorter DNA reads
    model: machine learning model developed by ONT (download from ONT bonito repository)
    """

    filepath = os.path.join(result_directory, outfile)
    batch = ''
    chunk = ''

    if batchsize > 0:
        batch = f'--batchsize {str(batchsize)}'
    if chunksize > 0:
        chunk = f'--chunksize {str(chunksize)}'

    bonito_call = f'bonito basecaller {batch} {chunk} {model} {directory} > {filepath}'
    command = f'conda run -n {env} {bonito_call}'
    subprocess.run(command, shell=True)

    return outfile


def chiron(directory, result_directory, folder='chiron', outfile='output.fastq', env='chiron'):
    """
    --DESCRIPTION--
    Same principle as bonito wrapper - activates conda environment before calling through subprocess
    Chiron will produce a folder containing items such as logs, results, and graphs

    --PARAMETERS--
    directory: directory path where all files to be basecalled are stored
    result_directory: result directory path - final basecalled sequence is stored here
    folder: name of the result folder created by Chiron
    outfile: name of the basecalled file
    env: name of the conda environment used to run Chiron
    """

    output_dir = os.path.join(result_directory, folder)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    chiron_call = f'chiron call -i {directory} -o {output_dir}'
    command = f'conda run -n {env} {chiron_call}'
    filepath = os.path.join(result_directory, outfile)
    subprocess.run(command, shell=True)

    chiron_result = natsorted(os.listdir(os.path.join(output_dir, 'result')))
    records = [SeqIO.read(os.path.join(output_dir, f'result/{file}'), 'fastq') for file in chiron_result]
    SeqIO.write(records, filepath, 'fastq')

    return outfile
