import importlib
import string

class Node:
    def __init__(self, package, function, path:string=None, variables:dict=None):
        self.package = package
        self.function = function
        self.path = path
        self.variables = variables

        if self.path == "" or self.path.isspace():
            self.path = None
        return

    def process(self, data):
        if self.path is None:
            data = getattr(importlib.import_module(self.package), self.function)(data, **self.variables)
        else:
            data = getattr(importlib.import_module(self.package, self.path), self.function)(data, **self.variables)
        return data
