#!/bin/sh

cwd=`pwd`

cd $LIBRADTRANDIR/share/libRadtran/examples

$LIBRADTRANDIR/bin/uvspec < UVSPEC_SIMPLE.INP 
#  310.000  5.416691e+01  3.971610e+01  1.145884e-14  4.310465e+00  6.159986e+00  2.179487e-15

cd $cwd

echo "End running libradtran"
