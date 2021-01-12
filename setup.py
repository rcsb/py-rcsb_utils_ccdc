# File: setup.py
# Date: 13-Dec-2020
#
# Updates:
#
#
import re

from setuptools import find_packages
from setuptools import setup

packages = []
thisPackage = "rcsb.utils.ccdc"

with open("rcsb/utils/ccdc/__init__.py", "r") as fd:
    version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', fd.read(), re.MULTILINE).group(1)

if not version:
    raise RuntimeError("Cannot find version information")

setup(
    name=thisPackage,
    version=version,
    description="RCSB Python utility classes support CCDC search",
    long_description="See:  README.md",
    author="John Westbrook",
    author_email="john.westbrook@rcsb.org",
    url="https://github.com/rcsb/py-rcsb_utils_ccdc.git",
    #
    license="Apache 2.0",
    classifiers=(
        "Development Status :: 3 - Alpha",
        # 'Development Status :: 5 - Production/Stable',
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ),
    entry_points={
        "console_scripts": [
            "ccdc_search_cli=rcsb.utils.ccdc.CcdcSearchExec:main",
        ]
    },
    #  The following is somewhat flakey --
    dependency_links=[],
    install_requires=[
        "mmcif >= 0.61",
        "rcsb.utils.io >= 0.95",
        "rcsb.utils.multiproc >= 0.18",
    ],
    packages=find_packages(exclude=["rcsb.mock-data", "rcsb.utils.tests-ccdc", "tests.*"]),
    package_data={
        # If any package contains *.md or *.rst ...  files, include them:
        "": ["*.md", "*.rst", "*.txt", "*.cfg"]
    },
    #
    test_suite="rcsb.utils.tests-ccdc",
    tests_require=["tox"],
    #
    # Not configured ...
    extras_require={"dev": ["check-manifest"], "test": ["coverage"]},
    # Added for
    command_options={"build_sphinx": {"project": ("setup.py", thisPackage), "version": ("setup.py", version), "release": ("setup.py", version)}},
    # This setting for namespace package support -
    zip_safe=False,
)
