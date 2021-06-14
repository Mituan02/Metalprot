# The software program is used to design Metal binding protein. 

# Installation
Create a python virtual environment before installing the package.
Activate the virtual environment and install the package by
```
conda activate env_conda

pip install -e .
```

# How to test a specific function in Metalprot.
pytest tests/test_core.py -k 'test_2ndshell'