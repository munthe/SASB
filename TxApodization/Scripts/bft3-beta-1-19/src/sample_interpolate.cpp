/**
 * @file   SampleInterpolate1D.cpp
 * @author Jens Munk Hansen <jmh@jmhpc.mt.elektro.dtu.dk>
 * @date   Fri Feb 12 11:15:06 2010
 * 
 * @brief  Sample Interpolation Module
 * 
 * $Id: sample_interpolate.cpp,v 1.5 2010-12-03 16:37:46 jmh Exp $
 *                                                             
 * $Author: jmh $                                              
 *                                                             
 * $Date: 2010-12-03 16:37:46 $                                
 *                                                             
 * $State: Exp $                                               
 *
 * $Name: beta-1-19 $
 *
 * $Log: sample_interpolate.cpp,v $
 * Revision 1.5  2010-12-03 16:37:46  jmh
 * Single precision works again
 *
 * Revision 1.4  2010/05/02 20:08:13  jmh
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
 * Revision 1.8  2010/03/21 12:20:27  jmh
 * Removed crappy fs,f0, and c
 *
 * Revision 1.7  2010/03/19 17:52:26  jmh
 * AV when beamforming t1.m
 *
 * Revision 1.6  2010/03/18 11:50:05  jmh
 * What happened to SampleInterpolate.cpp
 *
 * Revision 1.5  2010/03/17 00:19:40  jmh
 * *** empty log message ***
 *
 * Revision 1.4  2010/03/14 13:48:54  jmh
 * Alpha version, TODO: CppUnit, vectorization
 *
 * Revision 1.3  2010/03/13 13:07:33  jmh
 * Added interpolator
 *
 */

#ifndef SAMPLE_INTERPOLATE_CPP
#define SAMPLE_INTERPOLATE_CPP

#include <cmath>

#include <stdexcept>
using std::runtime_error;

#include <cstdlib> // abs


#include <cstdio> // debugging only

#include "sample_interpolate.h"

#ifdef _WIN32
# include <malloc.h>
#else
# include "mm_malloc.h"
#endif

#ifndef _NO_FIR
template <class Domain, class Range>
ALIGN16_BEGIN const Range SampleInterpolate<Domain, Range>::FIR_COEFF[FIR_LENGTH] ALIGN16_END =
  {0.00061294, 0.00093452, 0.0014247, 0.0018319,
   0.0019667, 0.0016124, 0.00057095, -0.001276,
   -0.0039057, -0.0071022, -0.01043, -0.013244,
   -0.014749, -0.014097, -0.010519, -0.0034654,
   0.0072554, 0.021389, 0.038194, 0.056487,
   0.074748, 0.091302, 0.10452, 0.11306,
   0.11601, 0.11306, 0.10452, 0.091302,
   0.074748, 0.056487, 0.038194, 0.021389,
   0.0072554, -0.0034654, -0.010519, -0.014097,
   -0.014749, -0.013244, -0.01043, -0.0071022,
   -0.0039057, -0.001276, 0.00057095, 0.0016124,
   0.0019667, 0.0018319, 0.0014247, 0.00093452,
   0.00061294};
#endif

//! Constructor (default)
template <class Domain, class Range> 
SampleInterpolate<Domain, Range>::SampleInterpolate() {}

//! Constructor
/*!
  A more elaborate description of the constructor
*/
template <class Domain, class Range>
SampleInterpolate<Domain, Range>::
SampleInterpolate(const Range *y, const Domain size) throw(std::runtime_error) {
  setData(y,size);
}

//! Destructor
/*!
  A more elaborate description of the destructor
*/
template <class Domain, class Range>
SampleInterpolate<Domain, Range>::~SampleInterpolate() {
  // Delete the y2Values
  if (curMethod == spline) {
    if (y2Values) {
      _mm_free(y2Values);
      y2Values = NULL;
    }
  }
}

/** 
 * Initialize data for interpolation
 * 
 * @param y Range data
 * 
 * @return void
 */
template <class Domain, class Range>
void SampleInterpolate<Domain, Range>::
setData(const Range *y, const Domain size) throw(std::runtime_error) {

  nElements = size;

  // Set the default interpolation method
  if (nElements == 0) {
    throw(runtime_error("SampleInterpolate::setData"
                        " ordinate is of zero length"));
  }
  else if (nElements == 1)
    curMethod = nearestNeighbour;
  else
    curMethod = linear;

  yValues = y;
}

/** 
 * 
 * 
 * @param x_req x-value
 * @param n order
 * @param offset offset in y-values (data used are y[offset]....y[offset+n-1]
 * 
 * @return interpolated value
 */
template <class Domain, class Range>
Range SampleInterpolate<Domain, Range>::
polynomialInterpolation(const Range x_req, Domain n, Domain offset) const throw(runtime_error) {
  /*
    A private function for doing polynomial interpolation
    Based on Nevilles Algorithm (Numerical Recipies 2nd ed., Section 3.1)
    x is the point we want to estimate, n is the number of points to use
    in the interpolation, and offset controls which n points are used
    (normally the nearest points)
  */

  // TODO: Test speed using static arrays (fixed order)

  // copy the x, y data into the working arrays
  Range* c = (Range*) _mm_malloc(n*sizeof(Range),16);
  Range* d = (Range*) _mm_malloc(n*sizeof(Range),16);
  Range* x = (Range*) _mm_malloc(n*sizeof(Range),16);

  Domain i;
  for (i = 0; i < n; i++) {
    d[i] = c[i] = yValues[offset];
    x[i] = Range(offset);
    offset++;
  }
  // Now do the interpolation using the rather opaque algorithm
  Range w, y;
  y = c[0];
  const Range one = Range(1);
  for (i = 1; i < n; i++) {
    // Calculate new C's and D's for each interation
    for (Domain j = 0; j < n-i; j++) {
      w = (c[j+1] - d[j]) * (one / (x[j] - x[j+i]));
      c[j] = (x[j] - x_req) * w;
      d[j] = (x[j+i] - x_req) * w;
    }
    y += c[0];
  }

  _mm_free(c);_mm_free(d);_mm_free(x);

  return y;
}

template <class Domain, class Range>
void SampleInterpolate<Domain, Range>::
setMethod(int newMethod) throw(std::runtime_error) {

  // Are we are switching to spline interpolation from something else?
  if (newMethod == spline && curMethod != spline) {
    // Calculate the y2Values

    y2Values = (Range*) _mm_malloc(nElements*sizeof(Range),16);

    /*
      The y2Values are initialised here.  I need to calculate the second
      derivates of the interpolating curve at each x_value.  As described
      in Numerical Recipies 2nd Ed. Sec. 3.3, this is done by requiring
      that the first derivative is continuous at each data point. This
      leads to a set of equations that has a tridiagonal form that can be
      solved using an order(N) algorithm.
      
      The first part of this solution is to do the Gaussian elimination so
      that all the coefficients on the diagonal are one, and zero below the
      diagonal.  Because the system is tridiagonal the only non-zero
      coefficients are in the diagonal immediately above the main
      one. These values are stored in y2Values temporarily. The temporary
      storage t, is used to hold the right hand side.
    */
    Range* t = (Range*) _mm_malloc(nElements*sizeof(Range),16);
    t[0] = 0;
    y2Values[0] = Range(0);
    
    if (nElements > 1)
      y2Values[nElements-1] = Range(0);

    Range a, b, delta;
    const Range six = 6;
    const Range one = Range(1.0);
    Range r;
    Domain i;
    for (i = 1; i < nElements-1; i++) {
      a = one;
      b = 2*(2); // 2*(*xValues[i+1] - *xValues[i-1]);
      r = (yValues[i+1] - yValues[i]) -
        (yValues[i] - yValues[i-1]);
      delta = t[i-1];
      /*
        if (nearAbs(b, delta))
        throw(runtime_error("SampleInterpolate::setMethod"
        " trouble constructing second derivatives"));
      */
      delta = b - delta;
      t[i] = one/delta;
      y2Values[i] = Range( (one/delta)*(six*r - a * (y2Values[i-1]) ) );
    }
    /*
      The second part of the solution is to do the back-substitution to
      iteratively obtain the second derivatives.
    */
    for (i = nElements-2; i > 1; i--)
      y2Values[i] = y2Values[i] - t[i]*(y2Values[i+1]);
    
    _mm_free(t);
  }
  else if ((curMethod == spline) && (newMethod != spline)) {
    // Delete the y2Values
    if (y2Values) {
      _mm_free(y2Values);
      y2Values = NULL;
    }
  }
  
  curMethod = newMethod;
}

template <class Domain, class Range> const Range* SampleInterpolate<Domain, Range>::
getY() const {
  return yValues;
}

// TODO: Remove inline
template <typename Domain, typename Range>
#ifdef _WIN32
Range SampleInterpolate<Domain, Range>::
#else
Range SampleInterpolate<Domain, Range>::
#endif
eval(Range x) const throw(std::runtime_error) {

  // Greater or equal
  Domain index = (Domain) std::max((float_type)ceil(x),Range(0)); 

  index = std::min(index,(Domain)nElements);

  Range x1,x2;
  Range y1,y2;

  switch (curMethod){
  
  case nearestNeighbour: // This does nearest neighbour interpolation
    {
      if (index == nElements)
        return yValues[nElements-1];
      else if (index == 0)
        return yValues[0];
      else {
        Range forward = fabs(Range(index) - x);
        index--;
        Range backward = fabs(Range(index) - x);
        if (forward < backward)
          index++;
        return yValues[index];
      }
    }
  
  case linear: // Linear interpolation is the default
    {
      if (index == nElements)
        index--;
      if (index == 0)
        index++;
      
      x2 = Range(index); y2 = yValues[index];
      index--;
      x1 = Range(index); y1 = yValues[index];

      return y1 + ((x-x1)/(x2-x1)) * (y2-y1);
    }

  case cubic:
    /*
      Fit a cubic polynomial to the four nearest points
      It is relatively simple to change this to any order polynomial
    */
    {
      if (nElements < 4)
        throw(runtime_error("SampleInterpolate::operator()"
                            " to few elements to do cubic polynomial interpolation"));
      if (index > 1 && index < nElements - 1)
        index  = index - 2;
      else if (index <= 1)
        index = 0;
      else
        index = nElements - 4;
      return polynomialInterpolation(x, (Domain) 4, index);
    }

  case spline: // Natural cubic splines
    {
      if (nElements < 4)
        throw(runtime_error("SampleInterpolate::operator()"
                            " to few elements to do spline interpolation"));
      if (index == nElements)
        index--;
      else if (index == 0)
        index++;
      Range h, a, b;
      Range y1d, y2d;
      const Range one = Range(1.0);
      const Range six = Range(6.0);

      x2 = Range(index); y2 = yValues[index]; y2d = y2Values[index];
      index--;
      x1 = Range(index); y1 = yValues[index]; y1d = y2Values[index];
      a = (x2-x);
      b = one-a;
      h = static_cast<Range>(one/six);
      return a*y1 + b*y2 + h*(a*a*a-a)*y1d + h*(b*b*b-b)*y2d;
    }
#ifndef _NO_FIR
  case fir: // FIR filter with 8 times interpolation / upsampling
    {
      
      Range result;
      ALIGN16_BEGIN Range result_tmp[2] ALIGN16_END;
      
      result_tmp[0] = 0.0;
      result_tmp[1] = 0.0;

      if (index == nElements) {
        return yValues[nElements-1];
      }
      else if (index == 0) {
        return yValues[0];
      }

      const Range* coef = (Range*) FIR_COEFF;                  // Fixed coefficients

      const size_t upsample_factor = FIR_INTERP_FACTOR; // Fixed upsample factor
      const size_t filter_length   = FIR_LENGTH;        // Fixed length
      const size_t filter_center   = FIR_CENTER;        // Fixed center

      Range I     = x * upsample_factor;
      Domain Iint = (Domain) floor(I);

      for (size_t i = 0 ; i < 2 ; i++) {
        for (size_t j = 0 ; j < filter_length ; j++) {
          
          Range indTmp = (Range)(Iint-filter_center+j+i) / upsample_factor;
          
          if( (indTmp == floor(indTmp))  &&  (indTmp >= 0)  &&  (indTmp < nElements) ) {
            result_tmp[i] = result_tmp[i] + (yValues[(Domain)indTmp])*upsample_factor*coef[j];
          }
        }
      }
      
      result = (Range)(result_tmp[0]*(Iint-I+1) + result_tmp[1]*(I-Iint));

      return result;
    }
#endif
  default:
    throw runtime_error("SampleInterpolate::operator() unknown type");
  }
}

// Hack
//template class SampleInterpolate<uint16_t,double>;
//template class SampleInterpolate<uint32_t,double>;
#ifdef M_SINGLE_PRECISION
template class SampleInterpolate<size_t,float>;
#else
template class SampleInterpolate<size_t,double>;
#endif

#endif
