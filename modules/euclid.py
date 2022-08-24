import numpy as np
import os
from Euclid.produceEuclidEncoding import makeFSM, trackGClevel, completelyRandom
from Euclid import convolutional
from Euclid import mapping


def euclidEncode(data, filepath, assignMech, segment_length:int =0):
    mechanism = completelyRandom
    if assignMech == 'track':
        mechanism = trackGClevel
    elif assignMech == 'random':
        mechanism = completelyRandom

    workSpaceDictionary = np.load(filepath, allow_pickle=True).item()
    c = workSpaceDictionary['c']
    vsym = workSpaceDictionary['vsym']
    hsym = workSpaceDictionary['hsym']
    numberOfPossibleCandidatesCountMatrix, outputDictionary, outputFSM, verticalSymbols, horizontalSymbols = makeFSM(c, vsym, hsym, mechanism)
    triggerLength = len(horizontalSymbols[0])
    symbolSize = len(vsym[0])
    initialState = '0'*symbolSize
    euclidFSM = convolutional.makeEuclidFSM(verticalSymbols, horizontalSymbols, outputFSM, initialState)

    data = ''.join(data.splitlines())
    if (len(data) % triggerLength) != 0:
        padding = '0' * (triggerLength - (len(data) % triggerLength))
        data = data + padding

    encodedStream = convolutional.FSMdictionaryEncoder(data, euclidFSM)
    flatStream = ''
    for sublist in encodedStream:
        flatStream = flatStream + sublist
    # Now we map the binary flat stream back into bases (A,C,T,G)
    codedData = mapping.binaryStreamToBases(flatStream)

    if segment_length > 0:
        output = [codedData[i:i + segment_length] for i in range(0, len(codedData), segment_length)]
        return output

    return codedData
