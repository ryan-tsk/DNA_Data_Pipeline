from node import Node
from utils import read_textfile, write_textfile

data = read_textfile("data_test.txt")
nodelist = []

compress_dict = {"json_output": "rev_codebook.json"}
new_node = Node("huffman_compression", "compress", package="custom", variables=compress_dict)
nodelist.append(new_node)

ecc_encode_dict = {"rate": "1/3"}
new_node = Node("ecc", 'encode', package="", variables=ecc_encode_dict)
nodelist.append(new_node)

binary_to_dna_dict = {}
new_node = Node("mapping", 'nt_mapping', package="", variables=binary_to_dna_dict)
nodelist.append(new_node)

dna_to_binary_dict = {"binary_to_nt": False}
new_node = Node("mapping", "nt_mapping", package="", variables=dna_to_binary_dict)
nodelist.append(new_node)

ecc_decode_dict = {"rate":"1/3"}
new_node = Node("ecc", "decode", package="", variables=ecc_decode_dict)
nodelist.append(new_node)

decompress_dict = {"json_input": "rev_codebook.json"}
new_node = Node("custom.huffman_compression", "decompress", package="", variables=decompress_dict)
nodelist.append(new_node)

for node in nodelist:
    data = node.process(data)

write_textfile("result_data.txt", data)
