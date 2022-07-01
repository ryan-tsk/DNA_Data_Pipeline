import heapq
import json
from utils import binary_to_text, text_to_binary


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


def generate_codebook(data, output_path):
    forward_codebook = {}
    reverse_codebook = {}
    frequencies = get_frequencies(data)
    heap = create_heap(frequencies)
    traverse(heapq.heappop(heap), "", forward_codebook, reverse_codebook)
    with open(output_path, 'w') as output:
        rev_json = json.dumps(reverse_codebook)
        output.write(rev_json)
    return forward_codebook


def compress(data, json_output, output_path=None):
    forward_codebook = generate_codebook(data, json_output)
    encoded = text_to_binary(data=data, output_path=output_path, codebook=forward_codebook)
    return encoded


def decompress(data, json_input, output_path=None):
    with open(json_input, 'r') as input:
        reverse_codebook = json.load(input)
    decoded = binary_to_text(data=data, output_path=output_path, codebook=reverse_codebook)
    return decoded
