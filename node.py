import importlib
import subprocess
import yaml

from abc import ABC, abstractmethod
from utils import read_textfile


class Node(ABC):
    def __init__(self, name, output_path=None, variables: dict = None):
        self.name = name
        self.output_path = output_path
        self.variables = variables

        if variables is None:
            self.variables = {}

    @abstractmethod
    def process(self, data):
        pass


class NodeFunction(Node):
    def __init__(self, name, module, function, output_path=None, package=None, variables: dict = None):
        super().__init__(name=name, output_path=output_path, variables=variables)
        self.data_process = self.__init_process(module=module, package=package, function=function)

    @staticmethod
    def __init_process(module, package, function):
        if module is None or function is None:
            return None
        if package is None:
            return getattr(importlib.import_module(module), function)
        module = '.' + module
        return getattr(importlib.import_module(module, package), function)

    def process(self, input_stream=None):
        print(self.name + ' node is processing...')
        output = self.data_process(input_stream, **self.variables)
        if self.output_path is not None:
            with open(self.output_path, 'w') as output_file:
                output_file.write(output)
        return output


class NodeCLICommand(Node):
    def __init__(self, name, output_path, command, variables: dict = None):
        super().__init__(name=name, output_path=output_path, variables=variables)
        self.command = command

    def process(self, data=None):
        subprocess.run(self.command, shell=True)
        return self.output_path


class NodeFactory:
    @staticmethod
    def create_node(properties, name):
        nodetype = properties['nodetype']
        del properties['nodetype']
        if nodetype == 'function':
            return NodeFunction(name=name, **properties)
        if nodetype == 'command':
            return NodeCLICommand(name=name, **properties)


class Pipeline:
    def __init__(self, config_path, data_path):
        self.nodes = []
        self.data_path = data_path
        self.build(config_path)

    def build(self, config_path):
        with open(config_path, 'r') as configFile:
            config = yaml.load(configFile, Loader=yaml.FullLoader)

        factory = NodeFactory()
        for item in config:
            node = factory.create_node(config[item], item)
            self.nodes.append(node)

    def run(self):
        data = read_textfile(self.data_path)
        for node in self.nodes:
            data = node.process(data)

    def test(self):
        pass
