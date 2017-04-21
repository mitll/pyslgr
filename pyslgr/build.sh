#!/bin/sh

python conda_setup.py install
dest_dir=${SP_DIR}/${PKG_NAME}
echo Copying Extension Modules to $dest_dir ...
cp -rL $PKG_NAME $dest_dir
# filtersSrc=`pwd`/../../filters
# filtersDest=${SP_DIR}/$PKG_NAME/filters
# echo Copying filters to $filtersDest
# cp -r $filtersSrc $filtersDest
