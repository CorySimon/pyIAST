from setuptools import setup

setup(name='pyiast',
    version='2.0',
    description='Ideal Adsorbed Solution Theory',
    url='https://github.com/CorySimon/pyIAST',
    download_url='https://github.com/CorySimon/pyIAST/tarball/master',
    install_requires=['numpy', 'scipy', 'pandas', 'matplotlib'],
    keywords='chemistry isotherm iast',
    author='Cory M. Simon',
    author_email='CoryMSimon@gmail.com',
    license='MIT',
    package_dir={'':'src'},
    packages=[''],
    zip_safe=False)
