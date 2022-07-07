import importlib
import subprocess


class Node:
    def __init__(self, module, function, package=None, variables: dict = None):
        self.module = module
        self.function = function
        self.package = package
        self.variables = variables

        if self.package is not None:
            if self.package == '' or self.package.isspace():
                self.package = None
                return

            if self.module[0] != '.':
                self.module = '.' + self.module

    def process(self, data):
        if self.package is None:
            data = getattr(importlib.import_module(self.module), self.function)(data, **self.variables)
        else:
            data = getattr(importlib.import_module(self.module, self.package), self.function)(data, **self.variables)
        return data


class NodeCLI:
    def __init__(self, command, input_path, output_path, variables: dict = None):
        self.command = command
        self.input_path = input_path
        self.output_path = output_path
        self.variables = variables

    def __clean_path__(self):
        if self.input_path == "" or self.input_path.isspace():
            self.input_path = None
        if self.output_path == "" or self.output_path.isspace():
            self.output_path = None
        return

    def process(self):
        subprocess.run(self.command, shell=True)
