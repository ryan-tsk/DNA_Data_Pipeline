from node import NodeFactory
import yaml

# command = 'bonito basecaller --batchsize 32 --fastq --weights 1 dna_r9.4.1_e8.1_sup@v3.3 results/stage2/ > results/stage2/output.fastq'
# test = NodeCLICommand('bonito', 'results/stage2/output.fastq', command)
# test.process()

with open('archived/config/config_stage3.yaml', 'r') as configFile:
    config = yaml.load(configFile, Loader=yaml.FullLoader)

factory = NodeFactory()
node = NodeFactory.create_node(config['Alignment'], 'Alignment')
node.process('output.fastq')
