import os.path
from pathlib import Path

def make_dir(path):
    if (not Path(path).exists()):
        os.makedirs(path)
    while path.find("//") != -1:
        path = path.replace('//', '/')
    if path[-1] != '/':
        path += '/'
    return path