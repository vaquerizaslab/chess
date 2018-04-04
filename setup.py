import os
from setuptools import setup, find_packages, Command


__version__ = None
exec(open('chess/version.py').read())


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
    name='chess',
    version=__version__,
    description='Comparison of Hi-C Experiments using Structural Similarity.',
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
)
