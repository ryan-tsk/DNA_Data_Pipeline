def nt_mapping(data, input_path=None, output_path=None, binary_to_nt: bool = True, nt_len: int = 300):
    if input_path is not None:
        with open(input_path, 'r') as infile:
            data = infile.readlines()

    if binary_to_nt:
        output = [binary_to_bases(seq.strip()) for seq in data]
    else:
        output = []
        for seq in data:
            if len(seq) > nt_len:
                seq = seq[0:nt_len]
            output.append(seq)

    if output_path is not None:
        with open(output_path, 'w') as outfile:
            outfile.write('\n'.join(output))

    print(output)
    return output


def binary_to_bases(bin_seq):
    """
        Converts binary string sequence to DNA nucleotide sequence.
    """
    mapping = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    return ''.join(mapping[bin_seq[i*2:i*2+2]] for i in range(len(bin_seq)//2))


def bases_to_binary(dna_seq):
    """
        Converts DNA nucleotide sequence to binary string sequence.
    """
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return ''.join(mapping[nt] for nt in dna_seq)
