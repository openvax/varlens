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
        description=(
            "commandline manipulation of genomic variants and NGS reads"),
        long_description=open('README.rst').read(),
        url="https://github.com/openvax/varlens",
        author="Tim O'Donnell",
        author_email="timodonnell@gmail.com",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
            'console_scripts': [
                'varlens-allele-support = varlens.commands.allele_support:run',
                'varlens-variants = varlens.commands.variants:run',
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
            'cython>=0.21',
            'numpy',
            'intervaltree',
            'pysam>=0.13',
            'typechecks',
            'varcode',
            'pyfaidx',
            'mhctools',
            'topiary',
        ],
    )
