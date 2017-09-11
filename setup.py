from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

install_requires = [
    'biopython>=1.70',
]

desc = """For a set of reference transposons, collect instances from one or more genomes and write as multi-FASTA files."""

setup(name='reunite',
      version='1.0.0',
      description=desc,
      url='https://github.com/Adamtaranto/TE-Reunite',
      author='Adam Taranto',
      author_email='adam.taranto@anu.edu.au',
      license='MIT',
      packages=['reunite'],
      classifiers=pypi_classifiers,
      keywords=["Transposon","TE","repeat","transposon"],
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'reunite=reunite.run_cmd:main',
        ],
    },
    )
