#AUTOMAKE_OPTIONS = gnu
AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies
VPATH = @srcdir@

include $(top_srcdir)/include

EXTRA_DIST  = $(top_srcdir)/src/@bft3_aperture $(top_srcdir)/src/@bft3_apodization $(top_srcdir)/src/@bft3_line $(top_srcdir)/src/@bft3_image $(top_srcdir)/src/@bft3_sampled_image $(top_srcdir)/src/@bft3_im_geom $(top_srcdir)/src/@bft3_system $(top_srcdir)/src/@bft3_emissions $(top_srcdir)/src/bft3_va_arg.m $(top_srcdir)/src/bft3_update_struct.m $(top_srcdir)/src/bft3_conv_float.m $(top_srcdir)/src/bft3_warn.m $(top_srcdir)/src/bft3_caller_name.m $(top_srcdir)/src/bft3_listfunctions.m $(top_srcdir)/src/bft3_lines.m $(top_srcdir)/src/bft3_apodizations.m $(top_srcdir)/examples/psf_8804p.m $(top_srcdir)/examples/psf_sasb_8804.m $(top_srcdir)/doc/html/* $(top_srcdir)/src/bft3_mex.mexa64 $(top_srcdir)/src/bft3_mex.mexglx $(top_srcdir)/doc/latex/refman.pdf

dist-hook:
	rm -rf `find $(distdir)/src/@bft3_aperture -name CVS`
	rm -rf `find $(distdir)/src/@bft3_apodization -name CVS`
	rm -rf `find $(distdir)/src/@bft3_line -name CVS`
	rm -rf `find $(distdir)/src/@bft3_image -name CVS`
	rm -rf `find $(distdir)/src/@bft3_sampled_image -name CVS`
	rm -rf `find $(distdir)/src/@bft3_system -name CVS`
	rm -rf `find $(distdir)/include -name CVS`
SUBDIRS = src

# Removing source
# rm -rf `find $(distdir)/src -name \*.c`
# rm -rf `find $(distdir)/src -name \*.cpp`
# rm -rf `find $(distdir)/include -name \*`
# rm -rf `find $(distdir) -name Makefile\*`
# rm -f `find $(distdir)/src/TODO`
# rm -f `find $(distdir)/INSTALL`
# rm -f `find $(distdir)/NEWS`
# rm -f `find $(distdir)/aclocal.m4`
# rm -f `find $(distdir)/configure*`
