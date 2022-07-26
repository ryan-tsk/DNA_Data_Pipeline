"""
This module is a simple attempt at simulating PCR by simply creating a fasta file that is
heavily weighted with the ground truth - the purpose of doing so is to have a fasta file that is heavily weighted
with the ground truth and is easy to test MSA and consensus generation

write_ground_truth - writes the ground truth of each read to a fasta file
write_weighted_fasta - takes the reads from the basecaller and merge them with the ground truth

weight determines the number of time the ground truth read is replicated
Example: weight = 100 means that the fasta file will have 1 basecalled read and 100 ground truth reads
"""

import os
from Bio import SeqIO
from utils import write_textfile, convert_to_fasta


def create_ground_truth(seqs, result_directory, filename: str, prefix_id: str = 'TEST'):
    ground_truth = []
    filepath = os.path.join(result_directory, filename)
    for i, seq in enumerate(seqs):
        seq_id = f'>{prefix_id}_{i}'
        ground_truth.extend([seq_id, seq])

    write_textfile(filepath, ground_truth, writelines=True)

    return seqs


def create_weighted_fasta(bc_filename, result_directory, gt_filename: str, folder: str = 'weighted',
                          weight: int = 100, prefix_id: str = 'TEST'):
    bc_path = os.path.join(result_directory, bc_filename)
    gt_path = os.path.join(result_directory, gt_filename)
    directory = os.path.join(result_directory, folder)

    if not os.path.exists(directory):
        os.mkdir(directory)

    if bc_path.endswith('fastq'):
        bc_path = convert_to_fasta(bc_path)

    basecall = SeqIO.parse(bc_path, 'fasta')
    ground_truth = SeqIO.parse(gt_path, 'fasta')

    bc_seqs = [str(record.seq) for record in basecall]
    gt_seqs = [str(record.seq) for record in ground_truth]

    for i, bc_seq in enumerate(bc_seqs):
        output = []
        seq_id = f'{prefix_id}_{i}'
        output.extend([f'>{seq_id}_{0}', bc_seq])

        for j in range(weight):
            output.extend([f'>{seq_id}_{j+1}', gt_seqs[i]])

        filepath = os.path.join(directory, f'{seq_id}.fasta')
        write_textfile(filepath, output, writelines=True)

    return directory
