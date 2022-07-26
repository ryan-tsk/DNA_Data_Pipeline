from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
# from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo

from utils import read_textfile

import subprocess
import os

# SeqIO.convert("example.fastq", "fastq", "example.fasta", "fasta")


def align_sequences(directory: str, wrapper: str, filepath: str, variables: dict = None):
    consensus_list = []

    if variables is None:
        variables = {}

    folder = os.listdir(directory)  #testing to see if this prevents tmp_alignment.fasta from being aligned
    output_path = os.path.join(directory, 'tmp_alignment.fasta')
    for file in folder:
        input_path = os.path.join(directory, file)
        if input_path == output_path:
            continue

        #need to ignore tmp_alignment.fasta
        if file.endswith('.fasta'):
            if wrapper == 'clustalw':
                    cline = ClustalwCommandline('clustalw', infile=input_path, outfile=output_path, output='FASTA', **variables)
            if wrapper == 'muscle':
                    cline = MuscleCommandline(input=input_path, out=output_path,  **variables)

            print(f'Aligning {str(file)}...')
            subprocess.run(str(cline), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            alignment = AlignIO.read(output_path, 'fasta')
            summary = AlignInfo.SummaryInfo(alignment)
            consensus = summary.dumb_consensus()
            consensus_list.append(str(consensus))

    with open(filepath, 'w') as txt:
        txt.write('\n'.join(consensus_list))

    return consensus_list


def create_consensus(directory: str):
    pass


def convert_to_fasta(in_file):
    filename = os.path.splitext(in_file)[0]
    new_file = filename + '.fasta'
    SeqIO.convert(in_file, 'fastq', new_file, 'fasta')
    return new_file


def clustalw_wrapper(in_file, variables:dict=None):
    # in_file is a fasta file
    if in_file.endswith('.fastq'):
        in_file = convert_to_fasta(in_file)

    cline = ClustalwCommandline("clustalw", infile=in_file, **variables)
    return cline
    #subprocess.run(str(cline), shell=True)


def muscle_wrapper(input, output, variables):
    if input.endswith('.fastq'):
        input = convert_to_fasta(input)

    cline = MuscleCommandline(input=input, out=output, **variables)
    return cline



