"""
Huffman compression adapted from : https://github.com/bhrigu123/huffman-coding
The original code has been modified to accommodate the DNA Storage Pipeline Project
"""

import heapq, os

from modules.mapping import binary_to_text, text_to_binary
from utils import write_json, load_json


class HuffmanNode:
    """
    Huffman node object used to create Huffman tree
    """
    def __init__(self, char, freq: int, left=None, right=None):
        self.char = char
        self.freq = freq
        self.left = left
        self.right = right

    def __lt__(self, other):
        return self.freq < other.freq

    def __eq__(self, other):
        if other is isinstance(other, HuffmanNode):
            return self.freq == other.freq
        return False


def get_frequencies(data):
    """
    --DESCRIPTION--
    Gets the frequencies of each character in the data stream

    --PARAMETERS--
    data: input string
    """
    frequencies = {}
    for char in data:
        if char not in frequencies:
            frequencies[char] = 1
        else:
            frequencies[char] += 1

    return frequencies


def create_heap(frequencies):
    """
    --DESCRIPTION--
    Create a heap (Huffman tree) using the character frequencies

    --PARAMETERS-
    frequencies: dictionary data item containing the frequency of each character
    """
    heap = []
    for key in frequencies:
        node = HuffmanNode(key, frequencies[key])
        heapq.heappush(heap, node)

    while len(heap) > 1:
        left = heapq.heappop(heap)
        right = heapq.heappop(heap)
        newNode = HuffmanNode(None, left.freq + right.freq, left, right)
        heapq.heappush(heap, newNode)

    return heap


def traverse(root, code, forward_codebook, reverse_codebook):
    """
    --DESCRIPTION--
    Recursive function, traversing through the heap

    --PARAMETERS--
    root: current heap node
    code: binary code of the current node
    forward codebook: dictionary storing the forward mapping of character to binary
    reverse codebook: dictionary storing the reverse mapping of binary to character
    """
    if root is None:
        return
    if root.char is not None:
        forward_codebook[root.char] = code
        reverse_codebook[code] = root.char
        return

    traverse(root.left, code + "0", forward_codebook, reverse_codebook)
    traverse(root.right, code + "1", forward_codebook, reverse_codebook)


def generate_codebook(data, codebook_path):
    """
    --DESCRIPTION--
    Recursive function, traversing through the heap

    --PARAMETERS--
    root: current heap node
    code: binary code of the current node
    forward codebook: dictionary storing the forward mapping of character to binary
    reverse codebook: dictionary storing the reverse mapping of binary to character
    """
    forward_codebook = {}
    reverse_codebook = {}

    frequencies = get_frequencies(data)
    heap = create_heap(frequencies)
    traverse(heapq.heappop(heap), "", forward_codebook, reverse_codebook)
    write_json(reverse_codebook, codebook_path)

    return forward_codebook


def compress(data, result_directory, reverse_filename: str, forward_filename: str = None):
    """
    --DESCRIPTION--
    Compression function used to perform Huffman compression

    --PARAMETERS--
    data: incoming ASCII data
    result_directory: directory where reverse_filename is saved
    reverse_filename: save reverse codebook to be used for decompression
    forward_filename: same as reverse filename but for forward - only use this if saving for reference
    """
    forward_codebook = generate_codebook(data, os.path.join(result_directory, reverse_filename))
    encoded = text_to_binary(data=data, codebook=forward_codebook)

    if forward_filename is not None:
        write_json(forward_codebook, os.path.join(result_directory, forward_filename))

    return encoded


def decompress(data, result_directory, reverse_filename):
    """
    --DESCRIPTION--
    Decompression function

    --PARAMETERS--
    data: incoming binary data
    result_directory: directory where reverse_filename is saved
    reverse_filename: load reverse codebook to be used for decompression
    """
    reverse_codebook = load_json(os.path.join(result_directory, reverse_filename))
    decoded = binary_to_text(data=data, codebook=reverse_codebook)
    return decoded
