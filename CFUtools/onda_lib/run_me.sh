#!/bin/bash
make distclean
aclocal -I /usr/share/aclocal -I./config
autoheader
automake -a # --copy
autoconf
./configure --with-matlab=/usr/local/matlabR2010b_64bit
#./configure --with-matlab=/usr/local/matlabR2010b_64bit --no-create --no-recursion
make
