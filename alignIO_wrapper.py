from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
# from Bio.Emboss.Applications import NeedleCommandline, WaterCommandline

from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
import subprocess
import os

# SeqIO.convert("example.fastq", "fastq", "example.fasta", "fasta")


def align_sequences(directory: str, wrapper: str, filepath: str, variables: dict = None):
    consensus = []
    for file in os.listdir(directory):
        output_path = os.path.join(directory, 'tmp_alignment.fasta')
        if file.endswith('.fasta'):
            if wrapper == 'clustalw':
                cline = ClustalwCommandline('clustalw', infile=file, outfile=output_path, **variables)
            if wrapper == 'muscle':
                cline = MuscleCommandline(input=file, **variables)
            subprocess.run(str(cline), shell=True)

        alignment = AlignIO.read(output_path, 'fasta')
        summary = AlignInfo.SummaryInfo(alignment)
        consensus = summary.dumb_consensus()
        consensus.append(str(consensus))

    with open(filepath, 'w') as txt:
        txt.write('\n'.join(consensus))



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

    if variables is None:
        variables = {}

    cline = ClustalwCommandline("clustalw", infile=in_file, **variables)
    return cline
    #subprocess.run(str(cline), shell=True)


def muscle_wrapper(in_file, out_file, variables):
    if in_file.endswith('.fastq'):
        in_file = convert_to_fasta(in_file)

    cline = MuscleCommandline(input=in_file, output=out_file, **variables)
    subprocess.run(str(cline), shell=True)



