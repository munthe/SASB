#AUTOMAKE_OPTIONS = gnu
AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies
VPATH = @srcdir@

include $(top_srcdir)/include

EXTRA_DIST  =  $(top_srcdir)/src/*.m $(top_srcdir)/doc/*.pdf $(top_srcdir)/src/onda_lib.mexa64

dist-hook:
	find $(distdir) -name .svn -exec rm {} \;
SUBDIRS = src
