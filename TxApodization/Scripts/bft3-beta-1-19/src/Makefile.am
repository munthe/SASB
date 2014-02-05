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

if ENABLE_CPP
  CPP_SRC = bft3_mex # aperture_mex apodization_mex line_mex image_mex sampled_image_mex
  TESTS = ApodizationTest InterpolateTest
endif

bin_PROGRAMS = $(CPP_SRC)
check_PROGRAMS = $(TESTS)

if ENABLE_CPP
  # aperture_mex_SOURCES = aperture_mex.cpp aperture.cpp mexarg.c mm_malloc.cpp   # aperture_mex_CFLAGS  = -DNeed_mex_gateway $(MATLAB_FLAGS)
  # aperture_mex_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 			       
  # aperture_mex_LDADD   = $(MATLAB_LIB)

  # apodization_mex_SOURCES = apodization_mex.cpp aperture.cpp apodization.cpp mexarg.c mm_malloc.cpp
  # apodization_mex_CFLAGS  = -DNeed_mex_gateway$(MATLAB_FLAGS)	       
  # apodization_mex_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 
  # apodization_mex_LDADD   = $(MATLAB_LIB)           

  # line_mex_SOURCES = line_mex.cpp aperture.cpp apodization.cpp line.cpp mexarg.c mm_malloc.cpp
  # line_mex_CFLAGS  = -DNeed_mex_gateway $(MATLAB_FLAGS)	       
  # line_mex_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 
  # line_mex_LDADD   = $(MATLAB_LIB)           

  # image_mex_SOURCES = image_mex.cpp aperture.cpp apodization.cpp line.cpp image.cpp sample_interpolate.cpp mexarg.c mm_malloc.cpp
  # image_mex_CFLAGS  = -DNeed_mex_gateway $(MATLAB_FLAGS)	       
  # image_mex_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 
  # image_mex_LDADD   = $(MATLAB_LIB)

  # sampled_image_mex_SOURCES = sampled_image_mex.cpp aperture.cpp apodization.cpp line.cpp sampled_image.cpp sample_interpolate.cpp mexarg.c mm_malloc.cpp
  # sampled_image_mex_CFLAGS  = -DNeed_mex_gateway $(MATLAB_FLAGS)	       
  # sampled_image_mex_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 
  # sampled_image_mex_LDADD   = $(MATLAB_LIB)

  bft3_mex_SOURCES = bft3_mex.cpp aperture_mex.cpp apodization_mex.cpp line_mex.cpp image_mex.cpp sampled_image_mex.cpp aperture.cpp apodization.cpp line.cpp image.cpp sampled_image.cpp sample_interpolate.cpp mexarg.c mm_malloc.cpp 
  bft3_mex_CFLAGS  = $(MATLAB_FLAGS)	       
  bft3_mex_LDFLAGS = $(WRAPS) $(MATLAB_LINK) 
  bft3_mex_LDADD   = $(MATLAB_LIB)

  InterpolateTest_SOURCES = sample_interpolate_app.cpp sample_interpolate.cpp mm_malloc.cpp sample_interpolate_test.cpp matlab_engine.cpp
  InterpolateTest_CFLAGS  = $(MATLAB_FLAGS)	       
  InterpolateTest_LDFLAGS = $(WRAPS) -L/usr/lib -L$(MATLAB_DIR)/bin/$(MATLAB_ARCH) -Wl,--rpath -Wl,$(MATLAB_DIR)/bin/$(MATLAB_ARCH)
  InterpolateTest_LDADD   = $(MATLAB_LIB) -leng -lcppunit -ldl

  ApodizationTest_SOURCES = apodization_app.cpp mm_malloc.cpp apodization_test.cpp
  ApodizationTest_CFLAGS   = $(MATLAB_FLAGS)	       
#  ApodizationTest_CXXFLAGS = $(CXXFLAGS)
  ApodizationTest_LDFLAGS  = $(WRAPS) -L/usr/lib -L$(MATLAB_DIR)/bin/$(MATLAB_ARCH) -Wl,--rpath -Wl,$(MATLAB_DIR)/bin/$(MATLAB_ARCH)
  ApodizationTest_LDADD    = $(MATLAB_LIB) -lcppunit -ldl
endif
