from pathlib import Path

import nbformat


notebook_path = Path(__file__).with_name("example_performance.ipynb")
notebook = nbformat.read(notebook_path, as_version=4)

for cell in notebook.cells:
    cell.metadata.pop("execution", None)

nbformat.write(notebook, notebook_path)
