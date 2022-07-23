import os
from Bio import SeqIO

# This is used to simulate PCR simply by replicating ground truth
# This results in a weighted fasta file which is used for MSA and consensus generation


# def create_ground_truth(seqs, directory: str, weight: int, id_prefix: str = 'TEST'):
#     for i, seq in enumerate(seqs):
#         seq_id = f'>{id_prefix}_{i}'
#         ground_truth = []
#         for j in range(weight):
#             ground_truth.extend([seq_id, seq])
#
#         filepath = os.path.join(directory, f'{seq_id}.fasta')
#         with open(filepath, 'w') as fasta:
#             fasta.write('\n'.join(ground_truth))
#
#     return seqs


def create_ground_truth(seqs, filepath: str, prefix_id: str = 'TEST'):
    ground_truth = []
    for i, seq in enumerate(seqs):
        seq_id = f'>{prefix_id}_{i}'
        ground_truth.extend([seq_id, seq])

    with open(filepath, 'w') as fasta:
        fasta.write('\n'.join(ground_truth))

    return seqs


def create_weighted_fasta(bc_path, gt_path, directory: str, weight: int, prefix_id: str = 'TEST'):
    if bc_path.endswith('fastq'):
        bc_path = convert_to_fasta(bc_path)

    basecall = SeqIO.parse(bc_path, 'fasta')
    ground_truth = SeqIO.parse(gt_path, 'fasta')

    bc_seqs = [str(record.seq) for record in basecall]
    gt_seqs = [str(record.seq) for record in ground_truth]

    for i, bc_seq in enumerate(bc_seqs):
        output = []
        seq_id = f'{prefix_id}_{i}'
        output.extend([f'>{seq_id}', bc_seq])

        for j in range(weight):
            output.extend([seq_id, gt_seqs[i]])

        filepath = os.path.join(directory, f'{seq_id}.fasta')
        with open(filepath, 'w') as fasta:
            fasta.write('\n'.join(output))

    return directory


def convert_to_fasta(in_file):
    filename = os.path.splitext(in_file)[0]
    new_file = filename + '.fasta'
    SeqIO.convert(in_file, 'fastq', new_file, 'fasta')
    return new_file
