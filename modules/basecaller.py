#Basecaller wrappers for Bonito and Chiron
from Bio import SeqIO

import os
import subprocess
from natsort import natsorted


def bonito(directory, result_directory, outfile='output.fastq', batchsize=0, chunksize=0,
           model='dna_r9.4.1_e8.1_sup@v3.3'):
    filepath = os.path.join(result_directory, outfile)
    batch = ''
    chunk = ''

    if batchsize > 0:
        batch = f'--batchsize {str(batchsize)}'
    if chunksize > 0:
        chunk = f'--chunksize {str(chunksize)}'

    command = f'bonito basecaller' + batch + chunk + f'{model} {directory} > {filepath}'
    subprocess.run(command, shell=True)

    return outfile


def chiron(directory, result_directory, folder='chiron', outfile='output.fastq', env='chiron'):
    output_dir = os.path.join(result_directory, folder)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    chiron_call = f'chiron call -i {directory} -o {output_dir}'
    command = f'conda run -n {env} {chiron_call}'
    filepath = os.path.join(result_directory, outfile)
    subprocess.run(command, shell=True)

    chiron_result = natsorted(os.listdir(os.path.join(output_dir, 'result')))
    records = [SeqIO.read(file, 'fastq') for file in chiron_result]
    SeqIO.write(records, filepath, 'fastq')

    return outfile
