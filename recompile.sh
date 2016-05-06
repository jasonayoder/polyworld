#!/bin/bash

module load python/2.7.8
module load gsl/1.15
module load qt/gnu/4.8.6

module unload opengl

make

cd src/app
cat .qt.mak_patch >.qt.mak
make -f .qt.mak
cd ../..

make

cd src/tools/PwMoviePlayer
cat .qt.mak_patch >.qt.mak
make -f .qt.mak
cd ../../..

make

touch bin/proputil
make

touch bin/pmvutil
make


module load opengl

module unload python/2.7.8
module unload gsl/1.15
module unload qt/gnu/4.8.6
