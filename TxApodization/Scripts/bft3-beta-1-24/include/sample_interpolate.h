/**
 * @file   SampleInterpolate.h
 * @author Jens Munk Hansen <jmh@jmhpc.mt.elektro.dtu.dk>
 * @date   Fri Feb 12 13:42:54 2010
 * 
 * @brief  Sample Interpolation Module (STL-version)
 * 
 * $Id: sample_interpolate.h,v 1.6 2010-05-02 20:08:13 jmh Exp $
 *                                                             
 * $Author: jmh $                                              
 *                                                             
 * $Date: 2010-05-02 20:08:13 $                                
 *                                                             
 * $State: Exp $                                               
 *
 * $Log: sample_interpolate.h,v $
 * Revision 1.6  2010-05-02 20:08:13  jmh
 * *** empty log message ***
 *
 * Revision 1.5  2010/04/07 14:15:28  jmh
 * *** empty log message ***
 *
 * Revision 1.4  2010/04/06 12:11:51  jmh
 * *** empty log message ***
 *
 * Revision 1.3  2010/03/30 19:34:15  jmh
 * VS XMT and RCV - slow but it works
 *
 * Revision 1.2  2010/03/30 15:46:24  jmh
 * *** empty log message ***
 *
 * Revision 1.1  2010/03/28 21:17:05  jmh
 * *** empty log message ***
 *
 * Revision 1.3  2010/03/21 12:20:27  jmh
 * Removed crappy fs,f0, and c
 *
 * Revision 1.2  2010/03/19 17:52:26  jmh
 * AV when beamforming t1.m
 *
 * Revision 1.1  2010/03/13 13:06:31  jmh
 * *** empty log message ***
 *
 * Revision 1.3  2010/02/12 17:32:10  jmh
 * Scalar version ready
 *
 * Revision 1.2  2010/02/12 12:39:37  jmh
 * STL version finished
 *
 * Revision 1.1  2010/02/12 12:36:31  jmh
 * *** empty log message ***
 *
 */

#ifndef SAMPLE_INTERPOLATE_H
#define SAMPLE_INTERPOLATE_H

#if (defined(_MSC_VER) && defined(_WIN32))
# define STATIC_INLINE_BEGIN static __forceinline 
# define STATIC_INLINE_END
#elif (defined(__GNUC__))
# define STATIC_INLINE_BEGIN static inline
# define STATIC_INLINE_END __attribute__ ((always_inline))
#endif

#ifdef _MSC_VER
// Visual C++ does not implement checked exceptions, so C4290 is just informing
// you that other exceptions may still be thrown from some functions
# pragma warning( disable : 4290 )
#endif

#ifdef _MSC_VER
 #define _USE_MATH_DEFINES
#endif

#include <cmath>

#include <stdexcept>
using std::runtime_error;

#include "common.h"

// Domain always integer type, Range always floating type
template <class Domain, class Range> class SampleInterpolate {
public:

  enum Method {
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

  SampleInterpolate();

  SampleInterpolate(const Range *y, const Domain size) throw(std::runtime_error);
  
  void setData(const Range *y, const Domain size) throw(std::runtime_error);
  
  virtual ~SampleInterpolate();

  int getMethod() const {return curMethod;}

  /** 
   * Set the method for interpolation,
   * SampleInterpolate<Domain, Range>::nearestNeighbour,
   * SampleInterpolate<Domain, Range>::linear,
   * SampleInterpolate<Domain, Range>::cubic,
   * SampleInterpolate<Domain, Range>::spline, or
   * SampleInterpolate<Domain, Range>::fir.
   * 
   * @param newMethod 
   */
  void setMethod(int newMethod) throw(std::runtime_error);

  virtual inline Range operator()(const Range x) const throw(std::runtime_error) {
    return this->eval(x);
  }

#ifdef _WIN32
  virtual __forceinline Range eval(Range x) const throw(std::runtime_error);
#else
  virtual Range eval(Range x) const throw(std::runtime_error);
#endif

  const Range* getY() const;

protected:
  int curMethod;         // Interpolation method to use
  Domain nElements;      // How many elements in the data set
  const Range*  yValues; // The corresponding ordinate of the data set
  Range* y2Values;       // The numerical second derivates (only for splines)

private:
  Range polynomialInterpolation(const Range x_req, Domain n,
                                Domain offset) const throw(runtime_error);

#ifndef _NO_FIR
  static const size_t FIR_LENGTH = 49;
  static const size_t FIR_CENTER = 24;
  static const size_t FIR_INTERP_FACTOR = 8;
  static const ALIGN16_BEGIN Range  FIR_COEFF[FIR_LENGTH] ALIGN16_END;
#endif

};

#if defined(_MSC_VER)
 #include "../src/sample_interpolate.cpp"
#endif

#endif
