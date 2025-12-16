import io
import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

DESCRIPTION = 'Hi-Compass: Depth-aware deep learning framework for cell-type-specific chromatin interaction prediction from ATAC-seq'

try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

NAME = 'hicompass'
EMAIL = '2247143021@qq.com'
URL = 'https://github.com/EndeavourSyc/Hi-Compass/'
AUTHOR = 'Yuanchen Sun'
VERSION = '1.0.1'

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    license='MIT',
    description=DESCRIPTION,
    url=URL,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-image",
        "pyBigWig",
        "cooler",
        "pysam",
        "piq",
        "matplotlib",
    ],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "hicompass=hicompass.cli:main",
        ],
    }
)
