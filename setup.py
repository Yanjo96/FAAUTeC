import os
import glob
import unittest
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

def my_test_suite():
    test_loader = unittest.TestLoader()
    test_suite = test_loader.discover('tests', pattern='*_test.py')
    return test_suite

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='TreeTopology',
    version='0.1',
    author='Nils Jenke; Yannick Hartmaring',
    #author_email='',
    description='Calculate if the optimal topology tree, calculated by e.g. RAxML, has a similar likelihood'
    +' as a topology tree followed certain contrains.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/yanjo96/treetopology',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3.7',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    keywords='topology tree, maximum likelihood',
    license='BSD',
    entry_points={
        'console_scripts': [
            'treetopology = treetopology.CLIOps:start_treetopology'
        ],
    },
    packages=['TreeTopology'], # So that the subfolder 'TreeTopology' is read immediately.
    #packages = find_packages(),
    install_requires=['biopython'],
    scripts=glob.glob('scripts/*'),
    test_suite='setup.my_test_suite',
    include_package_data=True,
    zip_safe=False
)
