#!/usr/bin/env python3
"""Package metadata for the paper companion; large databases are managed separately."""
from pathlib import Path
from setuptools import find_packages, setup

setup(
    name="metalprot",
    version="0.1.0rc1",
    author="Lei Lu",
    author_email="lonelur@gmail.com",
    url="https://github.com/Mituan02/Metalprot",
    description="Metal Mer-based search for local metal-binding sites",
    long_description=Path(__file__).with_name("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.9",
    packages=find_packages(include=["metalprot", "metalprot.*"]),
    package_data={"metalprot": ["constants/*.txt", "constants/*.pdb", "constants/*.pkl"]},
    install_requires=[
        "numpy>=1.23,<1.24", "matplotlib>=3.6,<4", "ProDy==2.4.1",
        "scikit-learn>=1.2,<1.7", "numba>=0.57,<0.61",
        "scipy>=1.10,<1.14", "pandas>=1.5,<2.3",
    ],
    extras_require={"test": ["pytest>=8,<9"], "build": ["build", "twine"]},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
    ],
)
