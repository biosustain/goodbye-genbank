# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import codecs
from setuptools import setup, find_packages

setup(
    name='gbgb',
    version='0.0.0',
    packages=find_packages(exclude=['*tests*', 'gbgb.export*']),
    url='https://github.com/biosustain/goodbye-genbank',
    license='Apache',
    author='Lars Schöning',
    author_email='lays@biosustain.dtu.dk',
    description='Goodbye, GenBank is a package for use with Biopython that gives feature annotations from '
                'GenBank records a new and better life.',
    long_description=codecs.open('README.rst', encoding='utf-8').read(),
)
