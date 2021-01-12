#!/bin/bash

# TODO

LIB_DIR=$PWD/../lib
TEST_DIR=$PWD
TESTS='unit_test.py'
PYTHONPATH_=$LIB_DIR

echo Lib is directory $LIB_DIR, contents are [$(ls $LIB_DIR)]
echo Testing in directory $TEST_DIR
echo Tests are [$TESTS]
echo Extra python paths are $PYTHONPATH_

flake8 --select F $TESTS # select pyflakes errors only


PYTHONPATH=hi echo hi $PYTHONPATH
PYTHONPATH=$LIB_DIR:$PYTHONPATH echo $PYTHONPATH
PYTHONPATH=$LIB_DIR:$PYTHONPATH mypy $TESTS
#pytest $TEST_DIR/$TESTS
