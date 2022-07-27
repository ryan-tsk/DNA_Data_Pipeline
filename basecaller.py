#Basecaller wrappers for Bonito and Chiron
from Bio import SeqIO

import os
import subprocess


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


def chiron(directory, result_directory, outfile='output.fastq', env='chiron'):
    command = f'chiron call -i {directory} -o {result_directory}'
    env_call = f'conda activate {env}'
    filepath = os.path.join(result_directory, outfile)

    subprocess.run(f'bash -c {env_call}', shell=True)
    subprocess.Popen([filepath, command], shell=True)

    chiron_directory = f'{result_directory}/result'

    records = [SeqIO.read(file, 'fastq') for file in os.listdir(chiron_directory)]
    SeqIO.write(records, filepath, 'fastq')

    return outfile