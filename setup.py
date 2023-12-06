try:
	from setuptools import setup, Extension
except ImportError:
	from distutils.core import setup
	from distutils.extension import Extension

import sys, platform
from setuptools import find_packages
import numpy as np

sys.path.append('python')

extra_compile_args = ['-DHAVE_KALLOC']
include_dirs = [".", np.get_include()]

if platform.machine() in ["aarch64", "arm64"]:
	include_dirs.append("sse2neon/")
	extra_compile_args.extend(['-ftree-vectorize', '-DKSW_SSE2_ONLY', '-D__SSE2__'])
else:
	extra_compile_args.append('-msse4.1') # WARNING: ancient x86_64 CPUs don't have SSE4



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
	ext_modules = [Extension('vacmap_index',
		sources = ['index/vacmap_index.pyx', 'index/bseq.c', 'index/seed.c', 'index/sketch.c', 'index/index.c', 'index/options.c',
				   'index/ksw2_extd2_sse.c', 'index/ksw2_exts2_sse.c', 'index/ksw2_extz2_sse.c', 'index/ksw2_ll_sse.c',
				   'index/kalloc.c', 'index/kthread.c', 'index/map.c', 'index/misc.c', 'index/sdust.c', 'index/esterr.c', 'index/splitidx.c'],
		depends = ['index/minimap.h', 'index/bseq.h', 'index/kalloc.h', 'index/kdq.h', 'index/khash.h', 'index/kseq.h', 'index/ksort.h',
				   'index/ksw2.h', 'index/kthread.h', 'index/kvec.h', 'index/mmpriv.h', 'index/sdust.h',
				   'index/cmappy.h', 'index/cmappy.pxd'],
		extra_compile_args = extra_compile_args,
		include_dirs = include_dirs,
		libraries = ['z', 'm', 'pthread'])],
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
