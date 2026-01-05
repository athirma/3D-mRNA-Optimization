#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from setuptools import setup, find_packages

setup(
    name="mrna_3d_optimization",
    version="1.0.0",
    description="Three-dimensional mRNA optimization algorithm incorporating codon pair bias",
    author_email="gongcheng@mail.tsinghua.edu.cn",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "biopython>=1.79",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
    ],
    python_requires=">=3.8",
)

