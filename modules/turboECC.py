"""
TurboDNA encoding adapted from : https://github.com/jasminequah/dna_archival_storage
The original code has been modified to accommodate the DNA Storage Pipeline Project
"""

from turboDNA import convolutional
import numpy as np


def set_variables(rate):
    """
    --DESCRIPTION--
    Sets the metadata for FSM creation

    Parameters
    ----------
    rate: Code rate for TurboDNA encoding (1/3, 1/2, 2/3)
    """
    if rate == '1/3':
        # ONE THIRD
        states = [0, 1, 2, 3, 4, 5, 6, 7]
        triggers = [[0], [1]]
        nextStateTable = np.array([[0, 4], [0, 4], [1, 5], [1, 5], [2, 6], [2, 6], [3, 7], [3, 7]])
        outputTable = [[[0, 0, 0], [1, 0, 0]],
                       [[0, 1, 1], [1, 1, 1]],
                       [[0, 1, 0], [1, 1, 0]],
                       [[0, 0, 1], [1, 0, 1]],
                       [[0, 0, 1], [1, 0, 1]],
                       [[0, 1, 0], [1, 1, 0]],
                       [[0, 1, 1], [1, 1, 1]],
                       [[0, 0, 0], [1, 0, 0]]]
        symbolSize = 3
        initialState = 0
        # bin_max_len = 100 for 150 nt, 200 for 300 nt

    elif rate == '1/2':
        # HALF
        states = [0, 1, 2, 3]
        triggers = [[0], [1]]
        nextStateTable = np.array([[0, 1], [2, 3], [0, 1], [2, 3]])
        outputTable = [[[0, 0], [1, 1]],
                       [[0, 1], [1, 0]],
                       [[1, 1], [0, 0]],
                       [[1, 0], [0, 1]]]
        symbolSize = 2
        initialState = 0
        # bin_max_len = 150 for 150nt, else 300 for 300 nt

    else:
        states = [0, 1, 2, 3, 4, 5, 6, 7]
        triggers = [[0, 0], [0, 1], [1, 0], [1, 1]]
        nextStateTable = np.array(
            [[0, 1, 2, 3], [4, 5, 6, 7], [1, 0, 3, 2], [5, 4, 7, 6], [2, 3, 0, 1], [6, 7, 4, 5], [3, 2, 1, 0],
             [7, 6, 5, 4]])
        outputTable = [[[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]],
                       [[0, 0, 1], [1, 0, 1], [0, 1, 1], [1, 1, 1]],
                       [[1, 0, 0], [0, 0, 0], [1, 1, 0], [0, 1, 0]],
                       [[1, 0, 1], [0, 0, 1], [1, 1, 1], [0, 1, 1]],
                       [[0, 1, 0], [1, 1, 0], [0, 0, 0], [1, 0, 0]],
                       [[0, 1, 1], [1, 1, 1], [0, 0, 1], [1, 0, 1]],
                       [[1, 1, 0], [0, 1, 0], [1, 0, 0], [0, 0, 0]],
                       [[1, 1, 1], [0, 1, 1], [1, 0, 1], [0, 0, 1]]]
        symbolSize = 3
        initialState = 0
        # bin_max_len = 200 for 150nt, else 400 for 300nt
    fsm = convolutional.FSM(states, triggers, outputTable, nextStateTable, initialState)
    return fsm, symbolSize


def encode(data, rate):
    """
    --DESCRIPTION--
    Sets the metadata for FSM creation

    --Parameters--
    data: input binary data to be encoded
    rate: Code rate for TurboDNA encoding (1/3, 1/2, 2/3)
    """

    fsm, symbol_size = set_variables(rate)
    bin_streams = [[int(bin) for bin in bin_seq] for bin_seq in data.split('\n')]

    encoded_streams = []
    for bin_stream in bin_streams:
        fsm.presentState = 0  # TODO: make a function in FSM to reset
        fsm.presentStateCoordinate = 0
        encoded_streams.append(convolutional.FSMEncoder(bin_stream, fsm))

    encoded_flatstreams = []
    for encoded_stream in encoded_streams:
        encoded_flatstream = []
        for sublist in encoded_stream:
            for item in sublist:
                encoded_flatstream.append(str(item))
        encoded_flatstreams.append(''.join(encoded_flatstream))

    return encoded_flatstreams


def decode(data, rate):
    """
    --DESCRIPTION--
    Decodes DNA nucleotide sequences using Viterbi decoder - returns decoded ASCII data


    --Parameters--
    ----------
    data: input binary data
    rate: Code rate for TurboDNA encoding (1/3, 1/2, 2/3)
    """

    fsm, symbolSize = set_variables(rate)

    def myFanOutFunction(state, observation, time):
        return convolutional.genericFanOutFunction(fsm, state, observation, time, None)

    decoded_streams = []
    for i, read_encoded_bin in enumerate(data):
        print("Decoding %s of %s..." % (i, len(data)))
        read_encoded_bin = [int(bin) for bin in read_encoded_bin]

        while len(read_encoded_bin) % 3 != 0:
            read_encoded_bin.pop()

        decoded_streams.append(
            convolutional.viterbiDecoderWithFlagging(8, 0, myFanOutFunction, read_encoded_bin, symbolSize,
                                                     produceGraphics=False)[0][0].pathTriggers)

    decoded_flatstreams = []
    for decoded_stream in decoded_streams:
        decoded_flatstream = []
        for sublist in decoded_stream:
            for item in sublist:
                decoded_flatstream.append(str(item))
        decoded_flatstreams.append(''.join(decoded_flatstream))

    return decoded_flatstreams
