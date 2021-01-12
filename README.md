# py-rcsb_utils_ccdc

RCSB Python wrapper for CCDC search and geometrical analysis tools.

## Introduction

Wrapper utilities for performing similarity and substructure searches on the CCDC database and
performing geometrical analysis of individual structure entries using the CCDC Python API.
This module has internal dependencies on the CCDC Python API and CCDC/CSD database.
These dependencies require licenses and separate installation that is described
with the CCDC documentation.  The 2021 version of the CCDC API is supported in Python 3.7.

### Installation

Download the library source software from the project repository and set the
enviroment corresponding to the CCDC installation in the the script `ccdc-api-env.sh`:

```bash

git clone --recurse-submodules https://github.com/rcsb/py-rcsb_utils_ccdc.git

cd py-rcsb_utils_ccdc
# edit and set the enviroment in the following bash script ...
. ccdc-api-env.sh

```

Optionally, run test suite (currently Python versions 3.7.9) using
[setuptools](https://setuptools.readthedocs.io/en/latest/) or
[tox](http://tox.readthedocs.io/en/latest/example/platform.html):

```bash
  # edit and set the enviroment in the following bash script ...
  . ccdc-api-env.sh
  pip install -r requirements.txt
  python setup.py test

or simply run:

 # edit and set the enviroment in the following bash script ...
  . ccdc-api-env.sh
  tox
```

Installation is via the program [pip](https://pypi.python.org/pypi/pip).  To run tests
from the source tree, the package must be installed in editable mode (i.e. -e):

```bash
pip install -e .
```

A CLI is provided to simplify access to the wrapper functions accross Python versions
not directly supported by the CCDC Python API.

```bash
# edit and set the enviroment in the following bash script ...
. ccdc-api-env.sh
python CcdcSearchExec.py --help
   -or-
ccdc_search_cli --help

usage: CcdcSearchExec.py [-h] [--mol_list_path MOL_LIST_PATH]
                         [--result_path RESULT_PATH]
                         [--search_type SEARCH_TYPE]
                         [--start_record START_RECORD]
                         [--end_record END_RECORD] [--csdhome CSDHOME]
                         [--python_lib_path PYTHON_LIB_PATH]
                         [--python_version PYTHON_VERSION]

optional arguments:
  -h, --help            show this help message and exit
  --mol_list_path MOL_LIST_PATH
                        Molecule file list path
  --result_path RESULT_PATH
                        Molecule file list path
  --search_type SEARCH_TYPE
                        Search type (similarity|substructure)
  --start_record START_RECORD
                        Starting record
  --end_record END_RECORD
                        End record
  --csdhome CSDHOME     Path to the CSD release (path to CSD_202x)
  --python_lib_path PYTHON_LIB_PATH
                        Path to Python library
  --python_version PYTHON_VERSION
                        Python library version (default: 3.7)

```
