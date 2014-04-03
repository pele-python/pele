#!/bin/bash
gtest=gtest-1.7.0
if [ ! -f $gtest.zip ]; then
  wget https://googletest.googlecode.com/files/$gtest.zip
fi
unzip $gtest.zip
mv $gtest source/gtest
