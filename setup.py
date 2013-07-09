"""
mage_tool
---------

mage_tool is a library and an executable used to create oligonucleotides for
use in multiplex automated genome engineering (MAGE) experiments.

"""

from setuptools import setup

setup(
    name='mage_tool',
    version='0.1',
    url='bio.ntz.nu/mage',
    license='Proprietary',
    author='Mads Valdemar Anderson, Michael Schantz Klausen',
    author_email='mskl@cbs.dtu.dk',
    description='A tool to automatically and easily create oligonucleotides for MAGE',
    long_description=__doc__,
    packages=['mage_tool', 'mage_tool.operations'],
    include_package_data=True,
    zip_safe=False,
    platforms='any',
    install_requires=[
        'biopython>=1.61',
        'pdfrw>=0.1',
        'PyYAML>=3.10',
        'matplotlib>=1.2.1',
        'reportlab>=2.7',
    ],
    scripts=[
        "MODEST.py",
        "MODESTcsv2report.py",
    ]
)
