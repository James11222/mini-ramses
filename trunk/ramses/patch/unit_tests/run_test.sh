#! /bin/bash

cp  ../../bin/Makefile ./
patch -R Makefile < Makefile_diff 

cd ../../bin/
make clean 
make -f ../patch/unit_tests/Makefile 
cd ../patch/unit_tests/

../../bin/ramses_unit_tests3d param_file.nml 

rm Makefile
