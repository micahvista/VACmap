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
