import subprocess
import os

from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline, ClustalOmegaCommandline, PrankCommandline
from Bio import AlignIO
from Bio.Align import AlignInfo

from utils import write_textfile
from natsort import natsorted


def align_sequences(directory: str, result_directory: str, wrapper: str, filename: str, variables: dict = None):
    consensus_list = []

    if variables is None:
        variables = {}

    tmpfile = 'tmp_alignment.fasta'
    output_path = os.path.join(directory, tmpfile)
    for file in natsorted(os.listdir(directory)):
        input_path = os.path.join(directory, file)
        if input_path == output_path:
            continue

        if file.endswith('.fasta'):
            if wrapper == 'clustalw':
                cline = ClustalwCommandline('clustalw', infile=input_path, outfile=output_path,
                                            output='FASTA', **variables)
            elif wrapper == 'muscle':
                cline = MuscleCommandline(input=input_path, out=output_path,  **variables)
            elif wrapper == 'clustalo':
                cline = ClustalOmegaCommandline(infile=input_path, outfile=output_path, **variables)
            elif wrapper == 'prank':
                cline = PrankCommandline(d=input_path, o=output_path, **variables)
            else:
                raise ValueError("No wrapper provided")

            print(f'Aligning {str(file)}...')
            subprocess.run(str(cline), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

            prankname = f'{output_path}.best.fas'
            if os.path.exists(prankname):
                os.rename(prankname, output_path)

            alignment = AlignIO.read(output_path, 'fasta')
            summary = AlignInfo.SummaryInfo(alignment)
            consensus = summary.dumb_consensus()
            consensus_list.append(str(consensus))
            os.remove(output_path)

    write_textfile(os.path.join(result_directory, filename), consensus_list, writelines=True)

    return consensus_list
