import yaml
from node import Node
from utils import read_textfile, write_textfile


with open('config.yaml', 'r') as configFile:
    config = yaml.load(configFile, Loader=yaml.FullLoader)

Nodes = []
for item in config:
    package = config[item]['package']
    module = config[item]['module']
    function = config[item]['function']
    variables = config[item]['variables']
    if variables is None:
        variables = {}

    node = Node(module=module, function=function, package=package, variables=variables)
    Nodes.append(node)

data = read_textfile("data_test.txt")
for node in Nodes:
    data = node.process(data)
