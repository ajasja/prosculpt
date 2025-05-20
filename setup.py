from setuptools import setup, find_packages

setup(
    name='prosculpt',
    version='0.1.0',
    url='https://github.com/ajasja/prosculpt',
    author='Ajasja Ljubetič',
    author_email='ajasja.ljubetic@gmail.com',
    description='Package includes functions to protein strucutre analysis and data manipulation necessary',
    packages=find_packages(),
    scripts=['rfdiff_mpnn_af2_merged.py'],
    install_requires=[],
)
