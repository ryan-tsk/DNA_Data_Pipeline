import importlib
import subprocess
import yaml
import os

from abc import ABC, abstractmethod
from utils import read_textfile, write_textfile, cleanup
from inspect import signature


class Node(ABC):
    def __init__(self, name, variables: dict = None):
        self.name = name
        self.variables = variables

        if variables is None:
            self.variables = {}

    @abstractmethod
    def process(self, input_stream, result_directory):
        pass


class NodeFunction(Node):
    def __init__(self, name, module, function, package=None, variables: dict = None):
        super().__init__(name=name, variables=variables)
        self.data_process = self.__init_process(module=module, package=package, function=function)

    @staticmethod
    def __init_process(module, package, function):
        if module is None or function is None:
            return None
        if package is None:
            return getattr(importlib.import_module(module), function)
        module = '.' + module
        return getattr(importlib.import_module(module, package), function)

    def process(self, input_stream, result_directory):
        print(self.name + ' node is processing...')
        sig = signature(self.data_process)
        if 'result_directory' in sig.parameters:
            output = self.data_process(input_stream, result_directory, **self.variables)
            return output

        output = self.data_process(input_stream, **self.variables)
        return output


class NodeCLICommand(Node):
    def __init__(self, name, outfile, command, variables: dict = None):
        super().__init__(name=name, variables=variables)
        self.command = command
        self.outfile = outfile

    def process(self, data=None, result_directory=None):
        subprocess.run(self.command, shell=True)
        return os.path.join(result_directory, self.outfile)


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
    def __init__(self, config_path):
        self.environment_name = 'Environment'
        self.input_path_name = 'Input Path'
        self.result_dir_name = 'Directory (Result)'
        self.log_dir_name = 'Directory (Log)'

        with open(config_path, 'r') as configFile:
            config = yaml.load(configFile, Loader=yaml.FullLoader)

        self._init_env_variables(config[self.environment_name])
        self.nodes = []
        self.build(config)

    def build(self, config):
        factory = NodeFactory()
        for item in config:
            if item == self.environment_name:
                continue
            node = factory.create_node(config[item], item)
            self.nodes.append(node)

    def run(self):
        data = read_textfile(self.data_path)
        for i, node in enumerate(self.nodes):
            log = [node.name, 'Input', str(data), '\n']
            data = node.process(data, self.result_directory)

            log.extend(['Output', str(data)])
            filename = f'{i}_{node.name}.txt'
            write_textfile(os.path.join(self.log_directory, filename), log, writelines=True)

    def _init_env_variables(self, variables: dict, result_folder='results', log_folder='logs'):
        self.data_path = variables[self.input_path_name]
        self.result_directory = os.path.join(result_folder, variables[self.result_dir_name])
        self.log_directory = os.path.join(log_folder, variables[self.log_dir_name])

        if not os.path.exists(self.result_directory):
            os.mkdir(self.result_directory)

        if not os.path.exists(self.log_directory):
            os.mkdir(self.log_directory)

        cleanup(self.result_directory)
        cleanup(self.log_directory)
