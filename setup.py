from __future__ import absolute_import
from setuptools import setup

setup(
    name='pyiast',
    version='1.4.3',
    description='Ideal Adsorbed Solution Theory',
    url='https://github.com/CorySimon/pyIAST',
    download_url='https://github.com/CorySimon/pyIAST/tarball/master',
    install_requires=['numpy', 'scipy', 'pandas>=0.24.0', 'matplotlib'],
    extras_require={'pre-commit': [
        "pre-commit==1.14.4",
        "yapf==0.26.0",
    ]},
    keywords='chemistry isotherm iast',
    author='Cory M. Simon',
    author_email='CoryMSimon@gmail.com',
    license='MIT',
    packages=['pyiast'],
    zip_safe=False)
