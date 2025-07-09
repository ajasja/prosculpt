from setuptools import setup, find_packages

setup(
    name='prosculpt',
    version='1.0.0',
    url='https://github.com/ajasja/prosculpt',
    author='Nej Bizjak',
    author_email='nejbizjak@gmail.com',
    description='Package includes functions to protein strucutre analysis and data manipulation necessary',
    packages=find_packages(),
    scripts=['rfdiff_mpnn_af2_merged.py', 'prosculpt.py', 'scripts/scoring_rg_charge_sap.py'],
    install_requires=['pandas','omegaconf','hydra-core', 'pathlib', 'biopython','numpy'],
)
