from utils import read_textfile

import os


def readfiles(dirpath):
    for files in os.listdir(dirpath):
        data = read_textfile(files)
        print(data)

readfiles('main_results')