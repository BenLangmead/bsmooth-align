#!/bin/sh

NM=bsmooth-align-`cat VERSION`

rm -rf .pkg
mkdir -p .pkg/$NM
cd .pkg/$NM

# Copy binaries
cp -r ../../bin .
chmod a+x bin/*

# Copy libraries
cp -r ../../lib .

# Copy example
cp -r ../../example .

# Copy index scripts
cp -r ../../index_scripts .

# Copy test
cp -r ../../test .

# Copy docs
cp ../../README .
cp ../../NEWS .
cp ../../COPYING .
cp ../../VERSION .

svn co svn+ssh://langmead@ibissub00.umiacs.umd.edu/fs/szdevel/src/svnroot/merman/trunk merman_trunk
mkdir merman
cp merman_trunk/* merman/
rm -rf merman_trunk

find . -name '.svn' | xargs rm -rf

cd ..
tar cvf - $NM | gzip -c > $NM.tar.gz
mv $NM.tar.gz ..
cd ..
