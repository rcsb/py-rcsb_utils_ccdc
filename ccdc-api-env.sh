
#!/bin/bash
export CSDHOME=/Applications/CCDC/CSD_2021
PY37=$PYENV_ROOT/versions/3.7.9

export LD_LIBRARY_PATH=$PY37/lib:$PY37/lib/python3.7/site-packages/ccdc/_lib:$LD_LIBRARY_PATH
#
export DYLD_LIBRARY_PATH=$PY37/lib/python3.7/site-packages/ccdc/_lib
export DYLD_FRAMEWORK_PATH=$PY37/lib/python3.7/site-packages/ccdc/_lib
#
pip install matplotlib jinja2 numpy biopython
#