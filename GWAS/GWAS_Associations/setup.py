from setuptools import setup

setup(
    name='R24',
    version='0.0.2',
    packages=["LD", "GWAS", "LiftOver", "Overlap", "Coloc"],
    url='',
    license='',
    package_dir = {"GWAS": "R24/GWAS", "LD": "R24/LD", "LiftOver": "R24/LiftOver", "Overlap": "R24/Overlap", "Coloc": "R24/Coloc", "MESC": "R24/MESC"},
    author='jrocha',
    author_email='jrocha@lji.org',
    description=''
)
