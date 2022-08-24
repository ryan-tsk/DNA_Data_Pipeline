import importlib
import subprocess
import yaml
import os
import logging

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
        sig = signature(self.data_process)
        if 'result_directory' in sig.parameters:
            output = self.data_process(input_stream, result_directory, **self.variables)
            return output

        output = self.data_process(input_stream, **self.variables)
        return output


class NodeCLI(Node):
    def __init__(self, name, outfile, command, variables: dict = None):
        super().__init__(name=name, variables=variables)
        self.command = command
        self.outfile = outfile

    def process(self, data=None, result_directory=None):
        subprocess.run(self.command, shell=True)
        return self.outfile


class NodeFactory:
    @staticmethod
    def create_node(properties, name):
        nodetype = properties['nodetype']
        del properties['nodetype']
        if nodetype == 'function':
            return NodeFunction(name=name, **properties)
        if nodetype == 'command':
            return NodeCLI(name=name, **properties)


class Pipeline:
    def __init__(self, config_path, input_path=None, output_path='results.txt', stages=False):
        self.input_path = input_path
        self.output_path = output_path
        self.stages = stages

        self.logger = logging.getLogger('Pipeline')

        with open(config_path, 'r') as configFile:
            self.config = yaml.load(configFile, Loader=yaml.FullLoader)

        environment = self.config.pop('Environment')
        self._init_environment(environment)
        self.nodes = []
        self.build(self.config)

    def build(self, config):
        factory = NodeFactory()
        for item in config:
            node = factory.create_node(config[item], item)
            self.nodes.append(node)

    def run(self):
        logging.basicConfig(filename=os.path.join(self.log_directory, 'log.txt'), filemode='a',
                            format="%(asctime)s %(message)s", level=logging.DEBUG)
        data = read_textfile(self.input_path)
        for i, node in enumerate(self.nodes):
            pre_data = data

            try:
                logging.info(f'{node.name} node is running...')
                data = node.process(data, self.result_directory)
            except Exception as e:
                if hasattr(e, 'message'):
                    message = e.message
                else:
                    message = e
                logging.error(f'Node {node.name} has failed. Cause: {message}')
                raise

            if self.stages:
                stage_result = [node.name, 'Input:', str(pre_data), '\n', 'Output:', str(data)]
                suffix = str(node.name).replace(' ', '_')
                filename = f'{i}_{suffix}.txt'
                write_textfile(os.path.join(self.stage_directory, filename), stage_result, writelines=True)

        write_textfile(os.path.join(self.result_directory, self.output_path), data)

    def _init_environment(self, variables: dict):
        if self.input_path is None:
            self.input_path = variables['input']

        self.result_directory = os.path.join(variables.get('results', 'results'), variables['name'])
        self.log_directory = os.path.join(variables.get('logs', 'logs'), variables['name'])

        if not os.path.exists(self.result_directory):
            os.mkdir(self.result_directory)

        if not os.path.exists(self.log_directory):
            os.mkdir(self.log_directory)

        if self.stages:
            self.stage_directory = os.path.join(variables.get('stages', 'stages'), variables['name'])
            if not os.path.exists(self.stage_directory):
                os.mkdir(self.stage_directory)
            cleanup(self.stage_directory)

        cleanup(self.result_directory)
        cleanup(self.log_directory)
