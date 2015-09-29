# Copyright (c) 2014. Mount Sinai School of Medicine
#
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

from setuptools import setup

readme_filename = "README.md"
current_directory = os.path.dirname(__file__)
readme_path = os.path.join(current_directory, readme_filename)

readme = ""
try:
    with open(readme_path, 'r') as f:
        readme = f.read()
except Exception as e:
    print(e)
    print("Failed to open %s" % readme_path)

try:
    import pypandoc
    readme = pypandoc.convert(readme, to='rst', format='md')
except Exception as e:
    print(e)
    print("Failed to convert %s from Markdown to reStructuredText"
        % readme_filename)

if __name__ == '__main__':
    setup(
        name='varlens',
        packages=["varlens", "varlens.commands", "varlens.read_evidence"],
        version="0.0.1",
        description=(
            "tools for looking at somatic variation across multiple samples"),
        long_description=readme,
        url="https://github.com/timodonnell/varlens",
        author="Tim O'Donnell",
        license="http://www.apache.org/licenses/LICENSE-2.0.html",
        entry_points={
            'console_scripts': [
                'varlens-allele-support = varlens.commands.allele_support:run',
                'varlens-variants = varlens.commands.variants:run',
                'varlens-reads = varlens.commands.reads:run',
                'varlens-fpkm2gsea = varlens.commands.fpkm2gsea:run',
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
            'varcode',
            'matplotlib',
        ],
    )
