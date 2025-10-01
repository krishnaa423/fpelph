from setuptools import setup, find_packages

setup(
    name='elph',
    version='0.0.1',
    description='Package for computing electron phonon matrix elements.',
    long_description='Package for computing electron phonon matrix elements.',
    author='Krishnaa Vadivel',
    author_email='krishnaa.vadivel@yale.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy', 
        'scipy', 
        'mpi4py',
        'h5py', 
        'ase',
        'xmltodict', 
        'setuptools',
        'jmespath',
    ],
    entry_points={
        'console_scripts': [
            'elph=elph.scripts.elph:elph',
        ],
    },
)
