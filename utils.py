import json
import os
import shutil
from Bio import SeqIO


def write_textfile(path, data, writelines=False):
    with open(path, 'w') as textfile:
        if writelines:
            textfile.write('\n'.join(data))
            return
        textfile.write(data)


def read_textfile(path, readlines=False):
    with open(path, 'r') as textfile:
        if readlines:
            data = textfile.readlines()
        else:
            data = textfile.read()
    return data


def write_json(json_data, path):
    with open(path, 'w') as output:
        rev_json = json.dumps(json_data)
        output.write(rev_json)


def load_json(path):
    with open(path, 'r') as input:
        data = json.load(input)
    return data


def convert_to_fasta(in_file):
    filename = os.path.splitext(in_file)[0]
    new_file = filename + '.fasta'
    SeqIO.convert(in_file, 'fastq', new_file, 'fasta')
    return new_file


def cleanup(directory):
    print(f'Cleaning {directory}')
    for file in os.listdir(directory):
        path = os.path.join(file, directory)
        print(f"Deleting path: {path}")
        if os.path.isfile(path) or os.path.islink(path):
            os.unlink(path)
        elif os.path.isdir(file):
            shutil.rmtree(file)

#make common util functions for read, write, data conversion, etc