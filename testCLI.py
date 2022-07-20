from node import NodeCLICommand

command = 'bonito basecaller --batchsize 32 --fastq --weights 1 dna_r10.4_e8.1_sup@v3.4 test_results/stage2/ > test_results/stage2/output.fastq'
test = NodeCLICommand('bonito', 'test_results/stage2/output.fastq', command)
test.process()
