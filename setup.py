# -*- encoding: utf-8 -*-
import io
import os
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import relpath
from os.path import splitext

from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    return io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ).read()


setup(
    name="gloria_soria_ddRAD_2015",
    version="0.1.0",
    license="MIT",
    description="Code supporting the Gloria-Soria et al 2015 paper on tsetse population genomics using ddRAD seq.",
    long_description="%s\n%s" % (read("README.rst"), re.sub(":obj:`~?(.*?)`", r"``\1``", read("CHANGELOG.rst"))),
    author="Gus Dunn",
    author_email="wadunn83@gmail.com",
    url="https://github.com/CacconeLabYale/gloria_soria_ddRAD_2015",
    packages=find_packages("src"),
    package_dir={"": "src"},
    py_modules=[splitext(basename(path))[0] for path in glob("src/*.py")],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT",
        "Operating System :: Unix",
        "Operating System :: POSIX",
        "Operating System :: Microsoft :: Windows",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Utilities",
    ],
    keywords=[
        # eg: "keyword1", "keyword2", "keyword3",
    ],
    install_requires=[
        # eg: "aspectlib==1.1.1", "six>=1.7",
        'click',
        'ipython',
        'sh',
        'matplotlib',
        'munch',
        'seaborn',
        'ggplot',
        'numpy',
        'pandas',
        'PyYAML',
        'tables',
        'scipy',
        'statsmodels',
        'ipdb',
        'pymc',
        'pybedtools',
        'tabulate',
        'scikits.bootstrap',
    ],
    extras_require={
        # eg: "rst": ["docutils>=0.11"],
    },
    dependency_links = [
        'https://github.com/xguse/spartan/archive/20151001.3.tar.gz#egg=spartan-20151001.3',
    ],
    entry_points={
        "console_scripts": [
            "gs_ddRAD2015 = gs_ddRAD2015.scripts.main:cli",

        ]
    },
)
