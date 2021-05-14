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
    name='FAAUTec',
    version='0.2',
    author='Yannick Hartmaring',
    #author_email='',
    description='Calculate if the AU Test for several Implementations using IQTree and RAxML for tree inference'
    +' as a topology tree followed certain contrains.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/Yanjo96/FAAUTeC',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3.8',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ],
    keywords='topology tree, maximum likelihood',
    license='BSD',
    entry_points={
        'console_scripts': [
            'faautec = faautec.CLIOps:start_faautec'
        ],
    },
    packages=['FAAUTeC'], # So that the subfolder 'TreeTopology' is read immediately.
    #packages = find_packages(),
    install_requires=['biopython','dendropy','ete3'],
    scripts=glob.glob('scripts/*'),
    test_suite='setup.my_test_suite',
    include_package_data=True,
    zip_safe=False
)
