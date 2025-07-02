# VACmap - a long-read aligner for structural variation discovery
# Copyright (C) 2023 Hongyu Ding
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

import sys, platform
from setuptools import find_packages
import numpy as np

sys.path.append('python')




setup(
	name = 'VACmap',
	version = '1.0',
	url = 'https://',
	description = 'VACmap',
	long_description = 'VACmap',
	author = 'Hongyu Ding',
	author_email = '920268622@qq.com',
	license = 'MIT',
	keywords = 'sequence-alignment',
	packages = find_packages("src"),
	package_dir = {"": "src"},
	classifiers = [
		'Development Status :: 5 - Production/Stable',
		'License :: OSI Approved :: MIT License',
		'Operating System :: POSIX',
		'Programming Language :: C',
		'Programming Language :: Cython',
		'Programming Language :: Python :: 3',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics'],
	setup_requires=["cython"],
	scripts=['src/vacmap/vacmap'])
