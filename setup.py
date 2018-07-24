import os
from setuptools import setup, find_packages, Command


__version__ = None
exec(open('chess/version.py').read())


with open("README.md", "r") as fh:
    long_description = fh.read()


class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info ./htmlcov')


setup(
    name='chess-hic',
    version=__version__,
    description='Comparison of Hi-C Experiments using Structural Similarity.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/vaquerizaslab/chess',
    setup_requires=[
        'setuptools>=18.0'
    ],
    packages=find_packages(),
    install_requires=[
        'numpy>=1.8.0',
        'scikit-image',
        'future',
        'pathos',
        'pandas',
        'scipy',
    ],
    scripts=['bin/chess'],
    cmdclass={
        'clean': CleanCommand
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
)
