AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies
VPATH = @srcdir@

INCLUDES = -I@top_srcdir@/include -I$(srcdir) -I$(MATLAB_DIR)/extern/include

# WRAPS = -Wl,--wrap,malloc,--wrap,calloc,--wrap,free,--wrap,assert

# You can implement this in a module and pass on wrap options
# void * __wrap_malloc(size_t size)
# {
# void *ptr;
# ptr = mxMalloc(size);
# if ( ptr != NUL )
# mexMakeMemoryPersistent(ptr);
# return ptr;
# }

# TODO: Remove -shared from unit test

bindir=$(exec_prefix)/bin

if WINNT
  MEXVSRC=$(MATLAB_DIR)/extern/include/mexversion.rc
  MEXVOBJ=$(top_builddir)/src/mexversionrc.o
  $(MEXVOBJ): $(MEXVSRC)
	windres -omexversionrc.o -O coff $(MEXVSRC)
endif

CPP_SRC = onda_lib

bin_PROGRAMS = $(CPP_SRC)
check_PROGRAMS = $(TESTS)

if ENABLE_CPP
  onda_lib_SOURCES = onda_lib.cpp onda_positioner.cpp tcp.cpp
  onda_lib_CFLAGS  = $(MATLAB_FLAGS)	       
  onda_lib_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 
  onda_lib_LDADD   = $(MATLAB_LIB)
endif
