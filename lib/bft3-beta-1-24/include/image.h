/******************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Image Class                                           *
 *                                                                            *
 * $Id: image.h,v 1.48 2011-05-01 14:54:13 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-05-01 14:54:13 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * $Log: image.h,v $
 * Revision 1.48  2011-05-01 14:54:13  jmh
 * *** empty log message ***
 *
 * Revision 1.47  2011/05/01 14:50:11  jmh
 * *** empty log message ***
 *
 * Revision 1.45  2011/04/17 20:09:06  jmh
 * Image representation class
 *
 * Revision 1.40  2011/04/17 14:54:21  jmh
 * Joinable threads
 *
 * Revision 1.35  2010/06/25 08:34:44  jmh
 * Fixed possibility to set and get window_parameter
 *
 * Revision 1.34  2010/05/05 09:29:39  jmh
 * No shared interpolator
 *
 *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Introduce common allocator 
 ****************************************************************************/

#ifndef IMAGE_H
#define IMAGE_H

#include <cstdio>

#ifdef _WIN32
# include <malloc.h>
#else
# include "mm_malloc.h"
#endif

#ifdef __GNUC__
# ifndef HAVE_THREAD
#  define HAVE_THREAD
# endif
# ifndef HAVE_PTHREAD_H
#  define HAVE_PTHREAD_H
# endif
#endif

#if (defined(_MSC_VER) && defined(_WIN32))
# define ALIGN16_BEGIN __declspec(align(16))
# define ALIGN16_END 
#elif defined(__GNUC__)
# define ALIGN16_BEGIN
# define ALIGN16_END __attribute__((aligned(16)))
#endif

#ifdef __GNUC__
#define __forceinline
#endif

#include <vector>
using std::vector;

// Forward declarations
class Aperture;

#include "apodization.h"
#include "line.h"
#include "common.h"

#include "mex_thread.h"

// Forward declarations
class Apodization;
class Line;
template <class Domain, class Range> class SampleInterpolate;

class Image {

  template <class T> class dereferenced_less {
   public:
    bool operator() (const T &lhs, const T &rhs) const {return *lhs < *rhs;}
  };

 public:

  typedef struct thread_arg {
    size_t start_line;
    size_t end_line;
    float_type* image;
    const float_type* rf_data;
    const float_type* delay;
    size_t n_rf_samples;
    const uint32_t* trans_no;
    size_t thread_id;
    int cpu_id;
  } thread_arg_t;

  enum Type {
    nearestNeighbour     = 0x00,
    linear               = 0x01,
    cubic                = 0x02,
    spline               = 0x03,
#ifndef _NOFIR
    fir                  = 0x04,
    n_types              = 0x05
#else
    n_types              = 0x04
#endif
  };

/**
 * Static variables
 * 
 */
  static size_t nthreads;
#ifdef HAVE_PTHREAD_H
  static pthread_t threads[N_MAX_THREADS];
  static pthread_attr_t attr;
# ifdef HAVE_MQUEUE_H
  static bool threads_initialized;
  static mqd_t mqd_master;
  static mqd_t mqd_client;
# endif
#endif
  static thread_arg_t threadarg[N_MAX_THREADS];

  Image();
  Image(const Aperture& Th_t, const Aperture& Th_r,
	const vector<Apodization>& Ah_t, const vector<Apodization>& Ah_r,
	const vector<Line>& lines);

  Image(const Image &other);

  /** 
   * Assignment
   * 
   * @param other 
   * 
   * @return 
   */
  Image& operator=(const Image &other);

  /** 
   * Destructor
   * 
   */
  ~Image();

#if defined(HAVE_PTHREAD_H)
  void* beamform_threaded(void *ptarg);
# ifdef HAVE_MQUEUE_H
  void* thread_func(void *ptarg);
# endif
#else
  unsigned __stdcall beamform_threaded(void *ptarg);
#endif

  bool Beamform(float_type* image, const float_type* rf_data,
		const float_type* delay, const size_t n_rf_samples,
		const uint32_t* trans_no);

 private:
 public:

  /** 
   * Clear POSIX message queue
   * 
   * @param qname name of the queue
   * 
   * @return 
   */
  static int mq_clear(const char* qname);

  class DataRep {
  public:

    // TODO: Use this for destroying threads when last object is freed
    static size_t ref_cnt;

    /**
     * Constructor
     *
     * @param xmt_aperture
     * @param rcv_aperture
     * @param xmt_apo
     * @param rcv_apo
     * @param lines
     * @param interp
     */
    DataRep(const Aperture& xmt_aperture, const Aperture& rcv_aperture,
  	    const vector<Apodization>& xmt_apo,
  	    const vector<Apodization>& rcv_apo,
  	    const vector<Line>& lines);
    /**
     * Destructor
     *
     */
    ~DataRep();
    Aperture m_xmt_aperture;
    Aperture m_rcv_aperture;

    /// Transmit apodization objects, one for each line
    vector<Apodization> m_xmt_apo;
    /// Receive apodization objects, one for each line
    vector<Apodization> m_rcv_apo;
    /// Lines
    vector<Line> m_lines;
    /// Interpolation type
    Type m_interp_type;
    /// Reference counter
    size_t m_cnt;
  };
  
  DataRep* data;
};

#endif
