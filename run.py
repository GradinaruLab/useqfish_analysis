import subprocess

from argparse import ArgumentParser
import os
from glob import glob

from warnings import filterwarnings

filterwarnings("ignore")

my_parser = ArgumentParser(description="Run analysis on a group of images")

my_parser.add_argument(
    "Path",
    metavar="path",
    type=str,
    help="the path to the diretory containing the images",
)


# Parse arguments
args = my_parser.parse_args()
path = args.Path

folders = os.listdir(path)
folders = sorted(folders)
print(folders)
folderpaths = [os.path.join(path, folder) for folder in folders]

for folderpath in folderpaths:
    if os.path.isdir(folderpath):
        os.system("python -m get_expression_matrix " + folderpath)
