"""
Huffman compression adapted from : https://github.com/bhrigu123/huffman-coding
The original code has been modified for the purpose of building a DNA data storage pipeline
"""

import heapq, os

from mapping import binary_to_text, text_to_binary
from utils import write_json, load_json


class HuffmanNode:
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
    frequencies = {}
    for char in data:
        if char not in frequencies:
            frequencies[char] = 1
        else:
            frequencies[char] += 1

    return frequencies


def create_heap(frequencies):
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
    if root is None:
        return
    if root.char is not None:
        forward_codebook[root.char] = code
        reverse_codebook[code] = root.char
        return

    traverse(root.left, code + "0", forward_codebook, reverse_codebook)
    traverse(root.right, code + "1", forward_codebook, reverse_codebook)


def generate_codebook(data, codebook_path):
    forward_codebook = {}
    reverse_codebook = {}

    frequencies = get_frequencies(data)
    heap = create_heap(frequencies)
    traverse(heapq.heappop(heap), "", forward_codebook, reverse_codebook)
    write_json(reverse_codebook, codebook_path)

    return forward_codebook


def compress(data, result_directory, reverse_filename: str, forward_filename: str = None):
    forward_codebook = generate_codebook(data, os.path.join(result_directory, reverse_filename))
    encoded = text_to_binary(data=data, codebook=forward_codebook)

    if forward_filename is not None:
        write_json(forward_codebook, os.path.join(result_directory, forward_filename))

    return encoded


def decompress(data, result_directory, reverse_filename):
    reverse_codebook = load_json(os.path.join(result_directory, reverse_filename))
    decoded = binary_to_text(data=data, codebook=reverse_codebook)
    return decoded
