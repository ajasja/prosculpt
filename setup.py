from setuptools import setup, find_packages

setup(
    name="prosculpt",
    version="1.1.4",
    url="https://github.com/ajasja/prosculpt",
    author="Ajasja Ljubetič",
    author_email="ajasja.ljubetic@gmail.com",
    description="Package includes functions to protein strucutre analysis and data manipulation necessary",
    packages=find_packages(),
    scripts=["prosculpt_run.py", "prosculpt.py", "scripts/scoring_script.py"],
    install_requires=[
        "pandas",
        "omegaconf",
        "hydra-core",
        "pathlib",
        "biopython",
        "numpy",
        "scipy",
    ],
)
