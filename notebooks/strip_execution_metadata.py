import argparse
from pathlib import Path

import nbformat



parser = argparse.ArgumentParser(description="Remove volatile execution timing from notebooks.")
parser.add_argument("notebooks", nargs="*", type=Path)
args = parser.parse_args()

notebook_paths = args.notebooks or [Path(__file__).with_name("example_performance.ipynb")]

for notebook_path in notebook_paths:
    notebook = nbformat.read(notebook_path, as_version=4)
    for cell in notebook.cells:
        cell.metadata.pop("execution", None)

    nbformat.write(notebook, notebook_path)
