#ifndef BFMATH_H
#define BFMATH_H

// Remove all this math crap
#ifdef _MSC_VER
 #define _USE_MATH_DEFINES
#endif

#include "common.h"

#include <cmath>
#include <cfloat>
#include <algorithm>

inline float_type dotprod(const float_type* a, const float_type* b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline void crossprod(float_type* a, float_type* b, float_type* c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}

// Dist point to line
// d = |(x0-x1)x(x0-x2)| / |x2-x1|
// point is x0, x1 and x2 lie on the line

/** 
 * Compute the distance from point (x0) to line passing through (x1) and (x2)
 * 
 * @param x0 
 * @param x1 
 * @param x2 
 * 
 * @return 
 */
inline float_type dist_point_to_line(const float_type* x0, const float_type* x1,
				     const float_type* x2) {
  ALIGN16_BEGIN float_type t1[3] ALIGN16_END;
  ALIGN16_BEGIN float_type t2[3] ALIGN16_END;
  ALIGN16_BEGIN float_type t3[3] ALIGN16_END;
  t3[0] = 0.0;t3[1]=0.0;t3[2]=0.0;

  t1[0] = x0[0] - x1[0];
  t1[1] = x0[1] - x1[1];
  t1[2] = x0[2] - x1[2];

  t2[0] = x0[0] - x2[0];
  t2[1] = x0[1] - x2[1];
  t2[2] = x0[2] - x2[2];

  crossprod(t1,t2,t3);
  float_type numerator = fabs(sqrt(dotprod(t3,t3)));

  t1[0] = x2[0] - x1[0];
  t1[1] = x2[1] - x1[1];
  t1[2] = x2[2] - x1[2];
#ifdef M_SINGLE_PRECISION
  float_type denominator = std::max((float_type)fabs(sqrt(dotprod(t1,t1))),
				    float_type(FLT_EPSILON));
#else
  float_type denominator = std::max((float_type)fabs(sqrt(dotprod(t1,t1))),
				    float_type(DBL_EPSILON));
#endif

  return (numerator / denominator);
}

/** 
 * Compute basis vectors from Euler angles (euler). The
 * (index)'th basis vector is returned as (output)
 * 
 * @param output  
 * @param euler 
 * @param index 
 */
inline void basis_vector(float_type* output, const float_type* euler,
			 size_t index) {

  const float_type alpha = euler[0];
  const float_type beta  = euler[1];
  const float_type gamma = euler[2];

  const float_type sa = sin(alpha);
  const float_type ca = cos(alpha);
  const float_type sb = sin(beta);
  const float_type cb = cos(beta);
  const float_type sc = sin(gamma);
  const float_type cc = cos(gamma);

  switch (index) {
  case 0:
    output[0]= ca*cc-cb*sa*sc;
    output[1]= cc*sa + ca*cb*sc;//-cb*cc*sa-ca*sc;
    output[2]= sb*sc; //sa*sb;
    break;
  case 1:
    output[0]= -ca*sc - cb*cc*sa;//cc*sa+ca*cb*sc;
    output[1]= ca*cb*cc-sa*sc;
    output[2]= cc*sb;//-ca*sb;
    break;
  case 2:
    output[0]= sa*sb;//sb*sc;
    output[1]= -ca*sb;//cc*sb;
    output[2]= cb;
    break;
  }
}

/** 
 * Compute basis vectors from Euler angles (euler). The
 * (index)'th basis vector is returned as (output)
 * 
 * @param output
 * @param input   
 * @param euler 
 */
inline void euler_rot(float_type* output, const float_type* input,
		      const float_type* euler) {

  const float_type alpha = euler[0];
  const float_type beta  = euler[1];
  const float_type gamma = euler[2];

  const float_type sa = sin(alpha);
  const float_type ca = cos(alpha);
  const float_type sb = sin(beta);
  const float_type cb = cos(beta);
  const float_type sc = sin(gamma);
  const float_type cc = cos(gamma);

  output[0] = input[0]*(ca*cc-cb*sa*sc) +
    input[1]*(-cb*cc*sa-ca*sc) +
    input[2]*sa*sb;
  output[1] = 
    input[0]*(cc*sa+ca*cb*sc) +
    input[1]*(ca*cb*cc-sa*sc) +
    input[2]*(-ca*sb);
  output[2] =
    input[0]*sb*sc +
    input[1]*cc*sb +
    input[2]*cb;
}

#endif
