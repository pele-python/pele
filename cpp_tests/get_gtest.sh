#!/bin/bash
# this script downloads the gtest source code and puts it into the folder
# source/gtest/
# The script is called by the travis build script .travis.yml
gtest=gtest-1.7.0
if [ ! -f $gtest.zip ]; then
  wget https://googletest.googlecode.com/files/$gtest.zip
fi
unzip $gtest.zip
mv $gtest gtest
