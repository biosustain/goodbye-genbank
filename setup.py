# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import codecs
from setuptools import setup, find_packages

setup(
    name='gbgb',
    version='0.1.0',
    packages=find_packages(exclude=['*tests*', 'gbgb.export*']),
    url='https://github.com/biosustain/goodbye-genbank',
    license='Apache',
    author='Lars Schöning',
    author_email='lays@biosustain.dtu.dk',
    description='Goodbye, GenBank is a package for use with Biopython that gives feature annotations from '
                'GenBank records a new and better life.',
    long_description=codecs.open('README.rst', encoding='utf-8').read(),
    install_requires=[
        'biopython>=1.66',
        'six>=1.8.0'
    ],
    test_suite='nose.collector',
    tests_require=[
        'nose>=1.1.2',
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
