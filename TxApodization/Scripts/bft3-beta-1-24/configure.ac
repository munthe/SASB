dnl  configure.ac
dnl
dnl  NOTE:  Process this file with GNU autoconf to produce a configure script.

# AC_INIT(FULL-PACKAGE-NAME, VERSION, BUG-REPORT-ADDRESS)
AC_INIT([bft3],[@PACKAGE_VERSION@],[jmh@elektro.dtu.dk],[])
AC_CONFIG_SRCDIR([src])
AC_CONFIG_AUX_DIR(config)

#AM_INIT_AUTOMAKE(@PACKAGE_NAME@, @PACKAGE_VERSION@ )
AM_INIT_AUTOMAKE(bft3, @PACKAGE_VERSION@ )
AM_PATH_CPPUNIT(1.12.1)

AM_CONFIG_HEADER(include/config.h)
AM_MAINTAINER_MODE
MEX_INIT([MEX],[1],[1],[2])
#PACKAGE_RELEASE_DATE="March 21 2010"
#AC_SUBST(PACKAGE_RELEASE_DATE)

dnl 
dnl  Build programs, C compiler, F77 compiler, make, install, etc.
dnl

AC_PROG_CXX([g++])
AC_PROG_CC([gcc cc])
AC_PROG_MAKE_SET
AC_PROG_INSTALL
dnl AC_PROG_RANLIB

dnl Disable CPP mex
AC_MSG_CHECKING(if cpp mex build is enabled)
AC_ARG_ENABLE(cpp,
[  --disable-cpp            disable cpp mex build],
  [if test "$enableval" = "yes" ; then
    enable_cpp=yes
  else
    enable_cpp=no
  fi], 
  enable_cpp=yes)

AM_CONDITIONAL(ENABLE_CPP, test "$enable_cpp" = "yes" )

if test "$enable_cpp" = "yes"
then

dnl
dnl  Threading are required
dnl
AC_DEFINE([HAVE_THREAD], 1, [Define to 1 if you have threading])

dnl
dnl  Headers
dnl

AC_STDC_HEADERS
AC_CHECK_HEADERS([pthread.h])
AC_CHECK_HEADERS([signal.h])
AC_CHECK_HEADERS([stdint.h])
AC_CHECK_HEADERS([stddef.h])
AC_CHECK_HEADERS([mqueue.h])

dnl
dnl Functions
dnl

AC_CHECK_FUNCS(printf)
AC_CHECK_FUNCS(pthread_getthreadid_np)

fi

CONFDATO=`date "+%Y-%m-%d %H:%M:%S"`
AC_DEFINE_UNQUOTED(CONFDATO,"$CONFDATO",[Date for running ./configure])

linux="no"
winnt="no"
case $build_os in
*linux*)
    linux="yes"
if test "${host_cpu}" = "x86_64";
then
    EXEEXT=.mexa64;
else
    EXEEXT=.mexglx;
fi;;
*mingw*)
    winnt="yes"
    EXEEXT=.dll;;
*cygwin*)
    winnt="yes"
    EXEEXT=.dll;;
esac

AM_CONDITIONAL(LINUX, test "$linux" = "yes")
AM_CONDITIONAL(WINNT, test "$winnt" = "yes")

dnl
dnl Matlab Library
dnl

CHECK_MATLAB

if test "$have_matlab" = "no"
then
    AC_MSG_ERROR([Must specify matlab directory using --with-matlab])
fi

dnl 
dnl  Libraries
dnl
#AC_CHECK_LIB([m]) # already added in CHECK_MATLAB

CXXFLAGS="$CXXFLAGS -fPIC -DMATLAB_MEX_FILE -DmxPOINTER_CLASS=mxUINT64_CLASS -DmxIsPointer=mxIsScalarUInt64 -DmxPOINTER_TYPE=UINT64_T -DmxIsPointerArray=mxIsUint64 -ansi -Wall -Werror -D_DEVELOPMENT"

# TODO: Use -ansi -pedantic -Wall
# ansi is equivalent to -std=c++98  -std=c++03
# Mathworks matlab.h doesn't conform with -std=c++0x

AC_CONFIG_FILES([Makefile \
                 src/Makefile])
AC_OUTPUT
