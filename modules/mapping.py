"""
Original code is by a previous MEng student: https://github.com/jasminequah/dna_archival_storage
This code has been modified to fit into the DNA data storage pipeline

When using with turboDNA encoder/decoder, nt_len (nt_mapping) and bin_len (text_to_binary)
will vary according to the setup used (see below)

The following assumes 300nt per sequence (if using 150nt, halve the following values)
If rate (turboDNA) = 1/3, nt_len = 300, bin_len = 200
If rate = 1/2, nt_len = 300, bin_len = 300
If rate = 2/3, nt_len = 300, bin_len = 400
"""

from utils import read_textfile


def nt_mapping(data, binary_to_nt: bool = True, nt_len: int = 300):
    """
    --DESCRIPTION--
    Converts text to binary

    --PARAMETER--
    data: input ASCII/binary data
    binary_to_nt: set to True if converting from binary to nt, False otherwise
    nt_len: nucleotide length that must be set if using TurboDNA - refer to earlier coderates for more information
    """

    output = []
    for seq in data:
        if binary_to_nt:
            output.append(binary_to_bases(seq.strip()))
        else:
            if len(seq) > nt_len:
                seq = seq[0: nt_len]
            output.append(bases_to_binary(seq.strip()))

    return output


def binary_to_bases(bin_seq):
    """
    --DESCRIPTION--
    Converts binary string sequence to DNA nucleotide sequence

    --PARAMETER--
    bin_seq: binary sequence
    """
    mapping = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}
    return ''.join(mapping[bin_seq[i*2:i*2+2]] for i in range(len(bin_seq)//2))


def bases_to_binary(dna_seq):
    """
    --DESCRIPTION--
    Converts DNA nucleotide sequence to binary string sequence

    --PARAMETER--
    dna_seq: DNA sequence
    """
    mapping = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    return ''.join(mapping[nt] for nt in dna_seq)


def text_to_binary(data, text_type: str = 'utf-8', bin_len: int = 200, codebook: dict = None):
    """
    --DESCRIPTION--
    Converts text to binary

    --PARAMETER--
    data: input ASCII data
    text_type: used to specify the text type for byte conversion
    bin_len: maximum binary length - this needs to be set if using
    codebook: only used for Huffman compression
    """

    bin_sequence = ''

    if codebook is None:
        bin_sequence = ''.join('{0:08b}'.format(ch, 'b') for ch in bytearray(data, text_type))

    else:
        for char in data:
            bin_sequence += codebook[char]

        padding = 8 - len(bin_sequence) % 8
        for i in range(padding):
            bin_sequence += "0"

        padding_info = "{0:08b}".format(padding)
        bin_sequence += padding_info

    bin_sequence = '\n'.join([bin_sequence[i: bin_len + i] for i in range(0, len(bin_sequence), bin_len)])

    return bin_sequence


def binary_to_text(data, codebook: dict = None):
    """
    --DESCRIPTION--
    Converts binary to text

    --PARAMETER--
    data: input binary data (string)
    input path: additionally supplementary input used to replace data if need to read from file
    codebook: only used for Huffman compression
    """

    text = ""
    data = "".join(data)

    if codebook is None:
        for i in range(0, len(data), 8):
            end = i + 8
            char = chr(int(data[i: end], 2))
            text += char

    else:
        padding = int(data[-8:], 2)
        last_char = - 1 * padding - 8
        encoded = data[:last_char]

        current_code = ""
        for bit in encoded:
            current_code += bit
            if current_code in codebook:
                char = codebook[current_code]
                text += char
                current_code = ""

    return text
