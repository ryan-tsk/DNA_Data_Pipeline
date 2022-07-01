import string


def text_to_binary(data, input_path=None, output_path=None, text_type: string = 'utf-8', nt_len: int = 200,
                   codebook: dict = None):
    if input_path is not None:
        data = read_textfile(input_path)
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

    bin_sequence = '\n'.join([bin_sequence[0 + i: nt_len + i] for i in range(0, len(bin_sequence), nt_len)])

    if output_path is not None:
        write_textfile(output_path, bin_sequence)
    return bin_sequence


def binary_to_text(data, input_path=None, output_path=None, codebook: dict = None):
    if input_path is not None:
        data = read_textfile(input_path)
        data = data.replace('\n', '')
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

    if output_path is not None:
        write_textfile(output_path, text)
    return text


def write_textfile(path, data):
    with open(path, 'w') as textfile:
        textfile.write(data)


def read_textfile(path):
    with open(path, 'r') as textfile:
        data = textfile.read()
    return data
