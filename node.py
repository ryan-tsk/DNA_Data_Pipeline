"""
Node - used to transform functions and executable commands into Node objects to be used with Pipeline
Pipeline - used to combine multiple Nodes to create an end-to-end pipeline that can be easily executed
"""


import importlib
import subprocess
import yaml
import os
import logging

from abc import ABC, abstractmethod
from utils import read_textfile, write_textfile, cleanup
from inspect import signature


class Node(ABC):
    """
    Parent class for NodeFunction and NodeCLI - New Node Classes can be extended from this class
    """

    def __init__(self, name, variables: dict = None):
        """
        --DESCRIPTION--
        Node initialisation

        --PARAMETERS--
        name: name of the Node
        variables: supplementary input variables (dictionary)
        """

        self.name = name
        self.variables = variables

        if variables is None:
            self.variables = {}

    @abstractmethod
    def process(self, input_stream, result_directory):
        """
        Abstract method used for inheritance
        """
        pass


class NodeFunction(Node):
    """
    Node object used to wrap around importable functions from other Python files
    NodeFunction inherits from Node
    """

    def __init__(self, name, module, function, package=None, variables: dict = None):
        """
        --DESCRIPTION--
        Initialise the NodeFunction object

        --PARAMETERS--
        name: name of the node
        module: name of Python file
        package: folder of Python file
        function: function name in Python file
        variables: supplementary input variables (optional)
        """

        super().__init__(name=name, variables=variables)
        self.data_process = self.__init_process(module=module, package=package, function=function)

    @staticmethod
    def __init_process(module, package, function):
        """
        --DESCRIPTION--
        Creates the process method by importing a function

        --PARAMETERS:
        module: name of Python file
        package: folder of Python file
        function: function name in Python file
        """
        if module is None or function is None:
            return None
        if package is None:
            return getattr(importlib.import_module(module), function)
        module = '.' + module
        return getattr(importlib.import_module(module, package), function)

    def process(self, input_stream, result_directory):
        """
        --DESCRIPTION--
        Callable process method - calling this method executes the imported process function

        --PARAMETERS--
        input_stream: input data stream
        result_directory: main result directory path
        """
        sig = signature(self.data_process)
        if 'result_directory' in sig.parameters:
            output = self.data_process(input_stream, result_directory, **self.variables)
            return output

        output = self.data_process(input_stream, **self.variables)
        return output


class NodeCLI(Node):
    """
    Node object used to wrap around executable CLI commands
    NodeCLI inherits from Node
    """

    def __init__(self, name, outfile, command, variables: dict = None):
        """
        --DESCRIPTION--
        Initialises the NodeCLI object

        --PARAMETERS--
        name: name of the node
        outfile: output file path
        command: command to be executed
        variables: supplementary variables (current not in use, for inheritance purposes)
        """

        super().__init__(name=name, variables=variables)
        self.command = command
        self.outfile = outfile

    def process(self, data=None, result_directory=None):
        """
        --DESCRIPTION--
        Callable process method - calling this method executes the command via subprocess
        Returns the outfile path

        --PARAMETERS--
        data: not in use, here for inheritance purposes (needs refactoring)
        result_directory: not in use, here for inheritance purposes  (needs refactoring)
        """
        subprocess.run(self.command, shell=True)
        return self.outfile


class NodeFactory:
    """
    Factory object used to easily create Node objects
    """

    @staticmethod
    def create_node(properties, name):
        """
        --DESCRIPTION--
        Factory method for creating Node objects

        --PARAMETERS--
        properties: the type of node that is created (only NodeFunction and NodeCLI right now)
        name: name of the Node
        """

        nodetype = properties['nodetype']
        del properties['nodetype']
        if nodetype == 'function':
            return NodeFunction(name=name, **properties)
        if nodetype == 'command':
            return NodeCLI(name=name, **properties)


class Pipeline:
    """
    Pipeline object used to run all nodes in a consecutive, end-to-end manner
    """

    def __init__(self, config_path, input_path=None, output_path='results.txt', stages=False):
        """
        --DESCRIPTION--
        Initialise the pipeline object

        --PARAMETERS--
        config_path: config file path
        input_path: input data file path (recommend to set in config, rather than in pipeline object)
        output_path: output file path
        stages: if set to True, all intermediate results between nodes are saved
        """

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
        """
        --DESCRIPTION--
        Builds the pipeline (list of Nodes) using a NodeFactory object

        --PARAMETERS--
        config: extracted configuration from config file
        """

        factory = NodeFactory()
        for item in config:
            node = factory.create_node(config[item], item)
            self.nodes.append(node)

    def run(self):
        """
        --DESCRIPTION:
        Runs each node consecutively in the order present in the configuration file
        Saves logs and stages (if enabled) after each execution
        """
        logging.basicConfig(filename=os.path.join(self.log_directory, 'log.txt'), filemode='a',
                            format="%(asctime)s %(message)s", level=logging.DEBUG)
        data = read_textfile(self.input_path)
        for i, node in enumerate(self.nodes):
            pre_data = data

            try:
                message = f'{node.name} node is running...'
                logging.info(message)
                print(message)
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
        """
        --DESCRIPTION--
        Used to initialise the pipeline environment, i.e. setting up log, result and stages directory

        --PARAMETERS--
        variables: extracted environment config from configuration file (dictionary)
        """

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
