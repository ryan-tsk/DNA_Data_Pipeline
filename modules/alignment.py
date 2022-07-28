import subprocess
import os

from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo

from utils import write_textfile
from natsort import natsorted


def align_sequences(directory: str, result_directory: str, wrapper: str, filename: str, variables: dict = None):
    consensus_list = []

    if variables is None:
        variables = {}

    output_path = os.path.join(directory, 'tmp_alignment.fasta')
    for file in natsorted(os.listdir(directory)):
        input_path = os.path.join(directory, file)
        if input_path == output_path:
            continue

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

    write_textfile(os.path.join(result_directory, filename), consensus_list, writelines=True)

    return consensus_list
