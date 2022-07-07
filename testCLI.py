from node import NodeCLI

command = 'bonito basecaller --batchsize 32 dna_r10.4_e8.1_sup@v3.4 bonito-test/ > bonito-test/output.fastq'
test = NodeCLI(command, input_path='', output_path='')
test.process()




