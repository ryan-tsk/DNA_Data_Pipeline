from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
# from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline

from Bio import SeqIO
from Bio import AlignIO
import subprocess
import os

# SeqIO.convert("example.fastq", "fastq", "example.fasta", "fasta")


def convert_to_fasta(in_file):
    filename = os.path.splitext(in_file)[0]
    new_file = filename + '.fasta'
    SeqIO.convert(in_file, 'fastq', new_file, 'fasta')
    return new_file


def clustalw_wrapper(in_file, variables:dict=None):
    # in_file is a fasta file
    if in_file.endswith('.fastq'):
        in_file = convert_to_fasta(in_file)

    if variables is None:
        variables = {}

    cline = ClustalwCommandline("clustalw", in_file=in_file, **variables)
    subprocess.run(str(cline), shell=True)


def muscle_wrapper(in_file, out_file, variables):
    if in_file.endswith('.fastq'):
        in_file = convert_to_fasta(in_file)

    cline = MuscleCommandline(input=in_file, output=out_file, **variables)
    subprocess.run(str(cline), shell=True)



