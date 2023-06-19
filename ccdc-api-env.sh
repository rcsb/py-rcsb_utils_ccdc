
#!/bin/bash
#
#  These settings are required to fire up the CCDC/CSD api...
#
#  Adjust these for settings for the site specific install of the
#  CSD database and the CCDC Python API
#
export CCDC=/data/wwpdb_da/da_top/reference/csd/cambridge_copy/
export CSDHOME=$CCDC/CSD_2022
export PYROOT=$CCDC/Python_API_2022/miniconda
export CSD_PYTHON_ROOT_PATH=$PYROOT

# For Linux
# export LD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
export LD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib:$PYROOT/lib:$LD_LIBRARY_PATH

# For MacOS
# export DYLD_LIBRARY_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
# export DYLD_FRAMEWORK_PATH=$PYROOT/lib/python3.7/site-packages/ccdc/_lib
#
# Packages recommended for some of the extra ccdc examples but not required for this package.
# pip install matplotlib jinja2 numpy biopython

# Activate miniconda package installed with full CSDS package to run python 3.7 conda env
source $PYROOT/bin/activate base

