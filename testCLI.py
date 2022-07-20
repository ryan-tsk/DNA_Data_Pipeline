from node import NodeCLICommand, NodeFunction, NodeFactory
import yaml

# command = 'bonito basecaller --batchsize 32 --fastq --weights 1 dna_r9.4.1_e8.1_sup@v3.3 test_results/stage2/ > test_results/stage2/output.fastq'
# test = NodeCLICommand('bonito', 'test_results/stage2/output.fastq', command)
# test.process()

with open('config/config_stage3.yaml', 'r') as configFile:
    config = yaml.load(configFile, Loader=yaml.FullLoader)

factory = NodeFactory()
node = NodeFactory.create_node(config['Alignment'], 'Alignment')
node.process('test_results/stage2/output.fastq')