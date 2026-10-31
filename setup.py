# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function
import os
import re

from setuptools import setup

current_directory = os.path.dirname(__file__)

with open('varlens/version.py', 'r') as f:
    version = re.search(
        r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
        f.read(),
        re.MULTILINE).group(1)

if __name__ == '__main__':
    setup(
        name='varlens',
        packages=["varlens", "varlens.commands", "varlens.read_evidence"],
        version=version,
        python_requires='>=3.11',
        description=(
            "commandline manipulation of genomic variants and NGS reads"),
        long_description=open('README.rst').read(),
        long_description_content_type='text/x-rst',
        url="https://github.com/openvax/varlens",
        author="Tim O'Donnell",
        author_email="timodonnell@gmail.com",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
            'console_scripts': [
                'varlens-allele-support = varlens.commands.allele_support:run',
                'varlens-variants = varlens.commands.variants:run',
                'varlens-fusions = varlens.commands.fusions:run',
                'varlens-reads = varlens.commands.reads:run',
            ],
        },
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Operating System :: OS Independent',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Programming Language :: Python',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
        install_requires=[
            'numpy>=2,<3',
            'pandas>=3,<4',
            'intervaltree>=3.1,<4',
            'pysam>=0.24.1,<0.25',
            'typechecks>=0.1.0',
            'varcode>=11.2.0,<12',
            'pyensembl>=2.27.2,<3',
            'pyfaidx>=0.9,<1',
            'mhctools>=4.1.1,<5',
            'topiary>=5.109.0,<6',
        ],
        extras_require={
            'fusion': ['isovar>=1.78.0,<2'],
            'test': ['pytest>=8,<10', 'pytest-cov>=5', 'ruff>=0.15'],
            'docs': ['sphinx>=8,<10', 'sphinx-rtd-theme>=3,<4',
                     'sphinxcontrib-autoprogram>=0.1.9',
                     'sphinxcontrib-programoutput>=0.20', 'numpydoc>=1.10'],
        },
    )
