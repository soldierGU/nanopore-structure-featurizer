from setuptools import setup, find_packages

setup(
    name="npstructfeat",
    version="0.1.0",
    description="Nanopore Structure Featurizer - Extract structural features from nanopore proteins",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.7",
    install_requires=[
        "numpy",
        "pandas",
        "biopython",
        "pyyaml",
        "matplotlib",
    ],
)