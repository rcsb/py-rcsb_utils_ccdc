
#!/bin/bash
#
#  These settings are required to fire up the CCDC/CSD api...
#
#  Adjust these for settings for the site specific install of the
#  CSD database and the CCDC Python API
#
export CSDHOME=/Applications/CCDC/CSD_2021
PYROOT=$PYENV_ROOT/versions/3.7.9
export CSD_PYTHON_PATH=$PYROOT/bin/python
#
# For Linux
export LD_LIBRARY_PATH=$PYROOT/lib:$PYROOT/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
# For Macos
export DYLD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
export DYLD_FRAMEWORK_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
#
#  Packages recommended for some of the extra ccdc examples but not required for this package.
# pip install matplotlib jinja2 numpy biopython
#