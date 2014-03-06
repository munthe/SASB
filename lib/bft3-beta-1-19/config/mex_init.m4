AC_DEFUN([MEX_INIT],
[
AC_PREREQ(2.54)
AC_CANONICAL_HOST
mex_var=`echo $1`_PLATFORM
AC_DEFINE_UNQUOTED($mex_var,"${host}")
mex_var=`echo $1`_MAJOR_VERSION
AC_DEFINE_UNQUOTED($mex_var,$2)
mex_var=`echo $1`_MINOR_VERSION
AC_DEFINE_UNQUOTED($mex_var,$3)
mex_var=`echo $1`_RELEASE_LEVEL
AC_DEFINE_UNQUOTED($mex_var,$4)

AC_ARG_WITH(debug,
[  --with-debug=yes	turn on debugging],
debug=$withval,
debug=no)

AC_MSG_CHECKING(debug mode)
if test $debug = yes
then
    AC_MSG_RESULT(yes)
else
    AC_MSG_RESULT(no)
fi

AC_ARG_WITH(docs,
[  --with-docs=yes      install doxumentation],
docs_install=$withval,
docs_install=no)

AC_MSG_CHECKING(install docs)
if test $docs_install = yes
then
    AC_MSG_RESULT(yes)
else
    AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(DOCS_INSTALL, test "x$docs_install" = "xyes")

AC_ARG_WITH(man,
[  --with-man=yes	install man page],
man_install=$withval,
man_install=yes)

AC_MSG_CHECKING(install man pages)
if test $man_install = yes
then
    AC_MSG_RESULT(yes)
else
    AC_MSG_RESULT(no)
fi
AM_CONDITIONAL(MAN_INSTALL, test "x$man_install" = "xyes")

AC_ARG_WITH(n32,
[  --with-n32=yes	use N32 ABI for SGI],
test "$withval" = "yes" && sgi_abi=n32)

AC_ARG_WITH(o32,
[  --with-o32=yes	use O32 ABI for SGI],
test "$withval" = "yes" && sgi_abi=o32)

AC_ARG_WITH(64,
[  --with-64=yes 	use 64 ABI for SGI],
test "$withval" = "yes" && sgi_abi=64)

AC_ARG_WITH(gnu,
[  --with-gnu=yes	use GNU versions of cc, f77],
gnu=$withval,
gnu=no)
])
