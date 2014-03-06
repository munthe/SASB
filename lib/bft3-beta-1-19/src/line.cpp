/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Line Class                                            *
 *                                                                            *
 * $Id: line.cpp,v 1.72 2011-07-27 22:27:34 jmh Exp $                         *
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-27 22:27:34 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * $Log: line.cpp,v $
 * Revision 1.72  2011-07-27 22:27:34  jmh
 * *** empty log message ***
 *
 * Revision 1.71  2011/07/19 20:17:25  jmh
 * Moved orientation to aperture class
 *
 * Revision 1.70  2011/07/17 21:49:10  jmh
 * JMH: Verify use of sampled_image!!!
 *
 * Revision 1.69  2011/05/30 09:17:50  jmh
 * *** empty log message ***
 *
 * Revision 1.68  2011-05-27 11:29:05  jmh
 * *** empty log message ***
 *
 * Revision 1.67  2011/05/01 14:54:13  jmh
 * *** empty log message ***
 *
 * Revision 1.63  2011/04/06 21:01:04  jmh
 * Added convex array
 *
 * Revision 1.61  2011/04/06 17:37:47  jmh
 * Divide by zero when apodize through ref
 *
 * Revision 1.57  2011-04-06 13:52:04  jmh
 * Error when using manual apodization
 *
 * Revision 1.51  2010/09/02 20:15:28  jmh
 * Removed line delays
 *
 * Revision 1.50  2010/09/02 19:42:58  jmh
 * Line delays
 *
 * Revision 1.48  2010/07/17 15:00:03  jmh
 * line apodization property
 *
 * Revision 1.44  2010/05/01 01:00:35  jmh
 * CTRL-C only works one time
 *
 * Revision 1.39  2010/04/09 14:15:39  jmh
 * Problem with fixed emit apodization
 *
 * Revision 1.37  2010/04/03 22:09:33  jmh
 * Fixed bug with VS
 *
 * Revision 1.35  2010/03/30 19:34:15  jmh
 * VS XMT and RCV - slow but it works
 *
 * Revision 1.33  2010/03/27 17:02:22  jmh
 * SEGV INTERPOLATE_MT ~/doc/article1/src/compound3sap.m
 *
 * Revision 1.32  2010/03/22 22:14:29  jmh
 * Apodization indices fixed, TODO: minsample
 *
 * Revision 1.31  2010/03/22 20:41:39  jmh
 * thread_template for WIN_32
 *
 * Revision 1.30  2010/03/21 12:20:27  jmh
 * Removed crappy fs,f0, and c
 *
 * Revision 1.29  2010/03/20 20:31:22  jmh
 * Dynamic and manual working
 *
 * Revision 1.28  2010/03/20 13:59:48  jmh
 * TODO: Verify why minsample isnt working
 *
 * Revision 1.27  2010/03/20 13:52:12  jmh
 * 2nd attempt dyn emit
 *
 * Revision 1.25  2010/03/19 19:57:09  jmh
 * First attempt dyn xmt apo
 *
 * Revision 1.24  2010/03/19 17:52:26  jmh
 * AV when beamforming t1.m
 *
 * Revision 1.21  2010/03/14 11:53:34  jmh
 * Pretty print
 *
 * Revision 1.20  2010/03/14 11:22:32  jmh
 * Doxygen comments
 *
 *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Find another way than pointer to a static variable
 *  - Use common headers
 *  - Find min sample if time (one way) = tstart
 *  - Use n_samples instead of length
 ****************************************************************************/

#ifdef HAVE_CONFIG_H
# define __STDC_LIMIT_MACROS 1
# include <stdint.h>
#endif

#include "aperture.h"
#include "apodization.h" // includes vector
#include "line.h"

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef _WIN32
# include <malloc.h>
#else
# include "mm_malloc.h"
#endif

#include <vector>
using std::vector;

#include <algorithm>
using std::sort;
using std::make_pair;
using std::pair;

#include "mexarg.h"

Line::DataRep::DataRep() : m_cnt(1), m_origin(NULL), m_direction(NULL) {
  m_dr     = float_type(0.0);
  m_length = float_type(0.0);
}

Line::DataRep::DataRep(const float_type* origin, const float_type* direction,
                       const float_type dr,
                       const float_type length) : m_cnt(1), m_dr(dr),
                                                  m_length(length)
{
  size_t i;
  m_origin    = (float_type*) _mm_malloc(4*sizeof(float_type),16);
  m_direction = (float_type*) _mm_malloc(4*sizeof(float_type),16);

  m_npos = (size_t) floor(m_length/m_dr);

  for (i=0;i<3;i++) {
    m_origin[i]    = origin[i];
    m_direction[i] = direction[i];
  }
  m_minsample = 0;
}

Line::DataRep::DataRep(const Apodization& xmt_apodization,
                       const Apodization& rcv_apodization,
                       const float_type* origin,
                       const float_type* direction, const float_type dr,
                       const float_type length) : m_cnt(1), m_dr(dr),
                                                  m_length(length) {
  size_t i;
  m_origin    = (float_type*) _mm_malloc(4*sizeof(float_type),16);
  m_direction = (float_type*) _mm_malloc(4*sizeof(float_type),16);

  m_npos = (size_t) floor(m_length/m_dr);

  for (i=0;i<3;i++) {
    m_origin[i]    = origin[i];
    m_direction[i] = direction[i];
  }
  m_xmt_apodization = xmt_apodization;
  m_rcv_apodization = rcv_apodization;

  m_minsample = 0;
                                                  }

Line::DataRep::~DataRep() {
  if (m_origin) {
    _mm_free(m_origin);
    m_origin = NULL;
  }
  if (m_direction) {
    _mm_free(m_direction);
    m_direction = NULL;
  }
}

Line::Line() : data(NULL) {
}

Line::Line(const Line &other) {
  data = other.data;
  if (data)
    data->m_cnt++;
}

Line& Line::operator=(const Line &other) {
  if (other.data)
    other.data->m_cnt++;

  if (data) {
    if (--data->m_cnt == 0) {
      delete data;
    }
  }

  data = other.data;
  return *this;
}

Line::Line(const float_type* origin, const float_type* direction,
           const float_type dr, const float_type length) {
  data = new Line::DataRep(origin, direction, dr, length);
}

Line::Line(const Apodization& xmt_apodization,
           const Apodization& rcv_apodization, const float_type* origin,
           const float_type* direction, const float_type dr,
           const float_type length) {
  data = new Line::DataRep(xmt_apodization,
                           rcv_apodization,
                           origin, direction, dr, length);
}

Line::~Line() {
  if (data) {
    if (--data->m_cnt == 0)
      delete data;
  }
}

bool Line::apodize(float_type* apo_values, int direction, size_t trans_no) {

	Aperture xmt_aperture;
	Aperture rcv_aperture;

	Apodization rcv_apodization = this->data->m_rcv_apodization;
	Apodization xmt_apodization = this->data->m_xmt_apodization;

  float_type pixel_x0, pixel_y0, pixel_z0;

  // Pixel location
  ALIGN16_BEGIN float_type pixel[3] ALIGN16_END;

  float_type dx, dy, dz;

	if ((!xmt_apodization.data) ||
			(!rcv_apodization.data)) {
		// One or more apodizations undefined
		return false;
	}
	else {
		xmt_aperture = xmt_apodization.data->m_aperture;
		rcv_aperture = rcv_apodization.data->m_aperture;
	}

  const size_t n_rcv_elements  = rcv_aperture.n_elements();

	findMinSample(float_type(0),xmt_aperture);

	// TODO: If desired you need to find distance from point on line to
	// nearest receive element
  this->data->m_minsample = 0;


	/*	mexPrintf("min. sample: %zu\n",data->m_minsample); */

	setApodization2(xmt_apodization,
									rcv_apodization, trans_no);
 
	pixel_x0 = origin(0);
	pixel_y0 = origin(1);
	pixel_z0 = origin(2);
	dx = this->dx();
	dy = this->dy();
	dz = this->dz();

	// Multiple focal points
	for (size_t k=0; (k+1) < data->m_limits.size() ; k++) {

		size_t line_start = data->m_limits[k];
		size_t line_end = data->m_limits[k+1];
		
		/*
		mexPrintf("line_start: %zu\n",line_start);
		mexPrintf("line_end: %zu\n",line_end);
		*/

		//		mexPrintf("xmt_apo_order: %zu, %zu\n",k, data->m_xmt_apo_order[k]);

		for (size_t j=line_start ; j < line_end ; j++) { // Samples along line
			
			pixel[0] = pixel_x0 + j*dx;
			pixel[1] = pixel_y0 + j*dy;
			pixel[2] = pixel_z0 + j*dz;
			
			if (rcv_apodization.dynamic()) {
				rcv_apodization.setPosition(pixel);
				rcv_apodization.setFnumber();
			}

			if (xmt_apodization.dynamic()) {
				xmt_apodization.setPosition(pixel);
				xmt_apodization.setFnumber();
			}
			
			float_type result;

			if (direction==0) {
				result = float_type(1.0);
				if (xmt_apodization.dynamic()) {
					result = result * xmt_apodization.dynApo(trans_no,pixel);
				}
				if (xmt_apodization.manual()) {
						result = result *
							xmt_apodization(trans_no,data->m_xmt_apo_order[k]);
				}
				apo_values[j] = result;
			}
			else {
				for (size_t l=0 ; l < n_rcv_elements ; l++) {
					
					result = float_type(1.0);
					
					if (rcv_apodization.dynamic()) {
						result = result * rcv_apodization.dynApo(l,pixel);;
					}
					if (rcv_apodization.manual()) {
						result = result *
							rcv_apodization(l,data->m_rcv_apo_order[k]);
					}
					apo_values[l+j*n_rcv_elements] = result;
				} // Receive elements
		  } // Receive apodization
		} // Samples along line
	} // Multiple focus points
	return true;
}

void Line::findMinSample(const float_type delay,
                         const Aperture& aperture) {

  // TODO: Will be static variable in the near future (single DLL)
  float_type c = *aperture.data->m_c;

  ALIGN16_BEGIN float_type v_lineref_to_aporef[3] ALIGN16_END;
  float_type dp_direction_rp;
  float_type dp_rp;
  float_type square_root_operand;

  // We need to find first pixel sample corresponding to time_offset
  v_lineref_to_aporef[0] = aperture.data->m_center_focus[0]-this->origin(0);
  v_lineref_to_aporef[1] = aperture.data->m_center_focus[1]-this->origin(1);
  v_lineref_to_aporef[2] = aperture.data->m_center_focus[2]-this->origin(2);

  dp_direction_rp = v_lineref_to_aporef[0]*this->direction(0) +
    v_lineref_to_aporef[1]*this->direction(1) +
    v_lineref_to_aporef[2]*this->direction(2);

  dp_rp = dotprod(v_lineref_to_aporef, v_lineref_to_aporef);
  
  /* Missing the case, where
  //                      -----------
  //                  ---/           \---
  //                -/                   \- t_offset
  //              -/                       \-
  //             /                           \
  //            /                             \
  //            |                             |       /----->
  //           /                        /-------------
  //           |          /-----o-------       |
  // origin /-------------                     /
  // -o-----    |           center_focus      |
  //            \                             /
  //             \                           /
  //              -\                       /-
  //                -\                   /-
  //                  ---\           /---
  //                      -----------
  // only data outside sphere, we cover the cases where we have a
  // positive and a negative solution, not to positive like the image
  */

  this->data->m_minsample = 0;
  
  // TODO: Check with /2
  float_type distance = delay*c/2;

  square_root_operand = SQUARE(dp_direction_rp) -       \
    (dp_rp - SQUARE(distance));
  if (!(square_root_operand > float_type(0.0))) {
    this->data->m_minsample = 0;
  }
  else {
    float_type dist_near = dp_direction_rp - sqrt(square_root_operand);
    float_type dist_far = dp_direction_rp + sqrt(square_root_operand);
    
    if ((dist_near < float_type(0.0)) && (dist_far > float_type(0.0))){
      this->data->m_minsample = (size_t) floor(dist_far / this->dr());
      // Data can never affect
      // point before this point
    }
  }
	// TODO: If desired you need to find distance from point on line to
	// nearest receive element
	this->data->m_minsample = 0;
}

// order means before (0) or after (1)

// BUG HERE
void Line::findVirtualPlane(const Apodization& apodization,
                            vector<size_t>& vp_limits,
                            vector<size_t>& order) {

  size_t max_index = (size_t) floor(this->length() / this->dr());

  const float_type* vs = apodization.data->m_aperture.data->m_focus;

  float_type md;

  float_type t_denom;
  float_type t;


  // This vector dotted with vector from virtual source to pixel
  // determines which sign to use (stupid)
  ALIGN16_BEGIN float_type normal[3] ALIGN16_END;
  normal[0] = apodization.data->m_aperture.data->m_center_focus[0] - 
    apodization.data->m_aperture.data->m_focus[0];
  normal[1] = apodization.data->m_aperture.data->m_center_focus[1] -
    apodization.data->m_aperture.data->m_focus[1];
  normal[2] = apodization.data->m_aperture.data->m_center_focus[2] -
    apodization.data->m_aperture.data->m_focus[2];

  vp_limits.push_back(0);

  md = dotprod(normal,vs);

  t_denom = dotprod(normal,this->data->m_direction);
  t = 0.0;

  t = (md - dotprod(this->data->m_origin,normal)) / t_denom;
  if (t < 0.0) {
    order.push_back(1); // after
  }
  else {
    order.push_back(0); // before
    // Compute index
    vp_limits.push_back(std::min( (size_t) ceil(t/this->dr()), max_index));
    order.push_back(1); // after
  }
  vp_limits.push_back(max_index);
}


// Works for xmt and rcv vs
void Line::setApodization2(const Apodization& xmt_apodization,
                           const Apodization& rcv_apodization,
                           const size_t trans_no) {

  size_t max_index = (size_t) floor(this->length() / this->dr());

  // Clear old sample limits and apodization indices
  data->m_xmt_apo_order.clear();
  data->m_rcv_apo_order.clear();
  data->m_xmt_vs_order.clear();
  data->m_rcv_vs_order.clear();

  data->m_limits.clear();

  // Temporary vectors for limits and apodization indices
  vector<size_t> xmt_order;
  vector<size_t> rcv_order;
  vector<size_t> xmt_limits;
  vector<size_t> rcv_limits;

  // Compute indices for apodization
  setAnyApodization(xmt_apodization, xmt_order, xmt_limits);
  setAnyApodization(rcv_apodization, rcv_order, rcv_limits);

  vector<size_t> xmt_vs_limits;
  vector<size_t> xmt_vs_order;

  vector<size_t> rcv_vs_limits;
  vector<size_t> rcv_vs_order;

  // Find virtual plane
  if (xmt_apodization.data->m_aperture.data->m_focus) {
    findVirtualPlane(xmt_apodization,xmt_vs_limits,xmt_vs_order);
  }
  else {
#if _DEVELOPMENT
    // TEST
    xmt_vs_limits.push_back(data->m_minsample);
    xmt_vs_limits.push_back(max_index);
    xmt_vs_order.push_back(1);    
#else
    xmt_vs_limits.push_back(1);
    xmt_vs_limits.push_back(max_index);
    xmt_vs_order.push_back(1);
#endif
  }

  // Find virtual plane
  if (rcv_apodization.data->m_aperture.data->m_focus) {
    findVirtualPlane(rcv_apodization,rcv_vs_limits,rcv_vs_order);
  }
  else {
#if _DEVELOPMENT
    rcv_vs_limits.push_back(data->m_minsample);
    rcv_vs_limits.push_back(max_index);
    rcv_vs_order.push_back(1);
#else
    rcv_vs_limits.push_back(1);
    rcv_vs_limits.push_back(max_index);
    rcv_vs_order.push_back(1);
#endif
  }

  // Sort indices and create sample apodization indices pairs
  vector<pair<size_t,int> > common_limits;

  // Apodization xmt
  std::transform(xmt_limits.begin(),xmt_limits.end(),
                 std::inserter(common_limits, common_limits.begin()),
                 sample_apodization_create(size_t(0)));

  // Apodization rcv
  std::transform(rcv_limits.begin(),rcv_limits.end(),
                 std::inserter(common_limits, common_limits.begin()),
                 sample_apodization_create(size_t(1)));

  // Virtual Source xmt  
  std::transform(xmt_vs_limits.begin(),xmt_vs_limits.end(),
                 std::inserter(common_limits, common_limits.begin()),
                 sample_apodization_create(size_t(2)));

  // Virtual Source rcv  
  std::transform(rcv_vs_limits.begin(),rcv_vs_limits.end(),
                 std::inserter(common_limits, common_limits.begin()),
                 sample_apodization_create(size_t(3)));

	/*

		// Work
  mexPrintf("Line limits for apodization\n");
  for (vector<size_t>::iterator itx= xmt_limits.begin();itx!=xmt_limits.end();itx++) {
    mexPrintf("XMT limit: %zu\n",*itx);
  }
	*/

  // Clear temporary vectors
  xmt_limits.clear();
  rcv_limits.clear();
  xmt_vs_limits.clear();
  rcv_vs_limits.clear();

  // Sorted by sample, then apodization in ascending order
  sort(common_limits.begin(),common_limits.end());

  vector<pair<size_t,int> >::iterator it;

#if 0
  mexPrintf("Line limits for apodization\n");
  for (it=common_limits.begin();it!=common_limits.end();it++) {
    mexPrintf("limit: %zu",it->first);
    if (it->second==0)
      mexPrintf(" XMT\n");
    else
      mexPrintf(" RCV\n");
  }
#endif

  vector<size_t>::iterator xmt_last_order = xmt_order.begin();
  vector<size_t>::iterator rcv_last_order = rcv_order.begin();
  vector<size_t>::iterator xmt_vs_last_order   = xmt_vs_order.begin();
  vector<size_t>::iterator rcv_vs_last_order   = rcv_vs_order.begin();

  it = common_limits.begin();

  size_t last_limit = SIZE_MAX;
  size_t next_limit = SIZE_MAX;

  bool init[4];
  init[0] = false;
  init[1] = false;
  init[2] = false;
  init[3] = false;

  vector<pair<size_t, int> >::iterator next;
  
  for (;it != common_limits.end() ; it++) {
    if (it->second == 0) {
      if (init[0])
        xmt_last_order++;
      init[0] = true;
    }
    else if (it->second == 1) {
      if (init[1])
        rcv_last_order++;
      init[1] = true;
    }
    else if (it->second == 2) {
      if (init[2])
        xmt_vs_last_order++;
      init[2] = true;
    }
    else if (it->second == 3) {
      if (init[3])
        rcv_vs_last_order++;
      init[3] = true;
    }

    next = it;
    next++;
    if (next!=common_limits.end())
      next_limit = next->first;
    else
      next_limit = SIZE_MAX;

    if ((it->first != last_limit) &&
        (it->first != next_limit)) {
      data->m_limits.push_back(it->first);
      last_limit = it->first;

      if (xmt_last_order != xmt_order.end())
        data->m_xmt_apo_order.push_back(*xmt_last_order);
      else
        data->m_xmt_apo_order.push_back(xmt_order.back());
      
      if (rcv_last_order != rcv_order.end())
        data->m_rcv_apo_order.push_back(*rcv_last_order);
      else
        data->m_rcv_apo_order.push_back(rcv_order.back());

      if (xmt_vs_last_order != xmt_vs_order.end())
        data->m_xmt_vs_order.push_back(*xmt_vs_last_order);
      else
        data->m_xmt_vs_order.push_back(xmt_vs_order.back());
      
      if (rcv_vs_last_order != rcv_vs_order.end())
        data->m_rcv_vs_order.push_back(*rcv_vs_last_order);
      else
        data->m_rcv_vs_order.push_back(rcv_vs_order.back());
    }
  }

#if 0
  vector<size_t>::iterator itc;
  vector<size_t>::iterator itc2;
  xmt_last_order = data->m_xmt_apo_order.begin();
  rcv_last_order = data->m_rcv_apo_order.begin();
  mexPrintf("Apodization indices:\n");
  for (itc = data->m_limits.begin(); itc!= data->m_limits.end();itc++) {
    itc2 = itc;
    itc2++;
    if (itc2==data->m_limits.end())
      continue;
    mexPrintf("[%zu ",*itc);
    if (itc2!=data->m_limits.end())
      mexPrintf("%zu] ",*itc2);
    if (xmt_last_order != data->m_xmt_apo_order.end())
      mexPrintf("XMT: %zu ",*xmt_last_order);
    if (rcv_last_order != data->m_rcv_apo_order.end())
      mexPrintf("RCV: %zu ",*rcv_last_order);
    xmt_last_order++;
    rcv_last_order++;
    mexPrintf("\n");
  }
#endif
  // Write down expression depth > F#*dyn_ap_size and solve for index on line

  /*
    if (xmt_apodization.dynamic()) {
    if (!xmt_apodization.data->m_aperture.data->m_focus) {
    ALIGN16_BEGIN float_type other_point_on_line[3] ALIGN16_END;
    other_point_on_line[0]=this->origin(0)+this->direction(0);
    other_point_on_line[1]=this->origin(1)+this->direction(1);
    other_point_on_line[2]=this->origin(2)+this->direction(2);
    float_type ortho_dist =
    dist_point_to_line(xmt_elem, xmt_apodization.data->m_ref, other_point_on_line)      
    }
    }
  */
}

// TODO: Fix no repeated range indices 
void Line::setAnyApodization(const Apodization& apodization,
                             vector<size_t>& order,
                             vector<size_t>& limits,
														 size_t i_emission) {
  
  // Vector from origin to apodization reference point
  ALIGN16_BEGIN float_type v_lineref_to_aporef[3] ALIGN16_END;
  ALIGN16_BEGIN float_type euler[3] ALIGN16_END;

  float_type dp_direction_rp;
  float_type dp_rp;
  float_type square_root_operand;

	/*
  const float_type c = *apodization.data->m_aperture.data->m_c;
  const float_type fs = *apodization.data->m_aperture.data->m_fs;
	*/

	/*
	mexPrintf("c: %f\n",c);
	mexPrintf("fs: %f\n",fs);
	*/

  // Max line index
  size_t max_index = (size_t) floor(this->length() / this->dr());

  // We need to find the intersection between each sphere belonging
  // to a distance and the current line. In the case, we have an
  // intersection, we need to find the index on the line
  // corresponding to that intersection. Since each intersection
  // gives at most 2 solutions, we divide the solutions into
  // "positive" and "negative" solutions and consequently sort the
  // intersections such that the right apodization can be used.

  /*       ^
  //       |                ^  (line.direction)
  //       |               /
  //       |         -----------
  //       |     ---/    /      \---
  //       |   -/       /           \-
  //       | -/        /              \-
  //       |/         X (line.origin)   \
  //       /           \                 \
  //       |            \ lineref_to_aporef
  //      /|             \                 \
  //      ||              o   (apodization ref)
  //      \|                               /
  //       +-------------------------->   |
  //      /\                             /
  //     /  \                           /
  //    /    -\                       /-
  //   /       -\                   /-
  //  /          ---\           /---
  // v               -----------
  //
  //  Diameter of circle equals first distance for apodizations
  */

  // The result is vector of limits of length #apodizations + 2 and
  // a vector of indices of length #apodizations + 1

  v_lineref_to_aporef[0] = apodization.ref(0)-this->origin(0);
  v_lineref_to_aporef[1] = apodization.ref(1)-this->origin(1);
  v_lineref_to_aporef[2] = apodization.ref(2)-this->origin(2);

  // We assume line direction is unit vector, hence if dx=dy=0, dz =
  // dr * 1, dr need not be unity
  dp_direction_rp = v_lineref_to_aporef[0]*this->direction(0) +
    v_lineref_to_aporef[1]*this->direction(1) + 
    v_lineref_to_aporef[2]*this->direction(2);

  dp_rp = dotprod(v_lineref_to_aporef, v_lineref_to_aporef);

  // STL vectors of distances
  vector<size_t> vec_near_distances;
  vector<size_t> vec_far_distances;

  // In the case some spheres do not intersect, TODO: sort those
  // spheres to find the one closest to the origin. This will be the
  // starting apodization.
  size_t skip_indices = 0;

  // TODO: Verify distances are positive
  float_type max_dist = float_type(0.0);
  size_t max_dist_index = 0;

  float_type dist_far, dist_near;
  size_t idist_far, idist_near;

  size_t m;
  for (m=0;m<apodization.data->m_ndistances;m++) {
    // Store max distance and corresponding index (in case of no
    // intersection)
    if (max_dist < apodization.distance(m)) {
      max_dist = apodization.distance(m);
      max_dist_index = m;
    }

    // Compute line-sphere intersection - two solutions at most
    square_root_operand = dp_direction_rp*dp_direction_rp - \
      (dp_rp - apodization.distance(m)*apodization.distance(m));
    if (square_root_operand < float_type(0.0)) {
      skip_indices++;
      continue;
    }
    dist_near = dp_direction_rp - sqrt(square_root_operand);
    dist_far = dp_direction_rp + sqrt(square_root_operand);

    // Compute corresponding line indices - consider rounding
    idist_near = (size_t)
      std::max( (float_type) floor(dist_near/this->dr()) , (float_type) 0.0);
    idist_far = (size_t)
      std::max( (float_type) floor(dist_far/this->dr()),(float_type)0.0);

    // One intersection only
    if (idist_near == idist_far) {
      vec_far_distances.push_back(idist_far);
      continue;
    }
    vec_near_distances.push_back(idist_near);
    vec_far_distances.push_back(idist_far);
  }

  vector<size_t*> vec_ref_near_distances; // Must cobe for sign
  vector<size_t*> vec_ref_far_distances;
  
  // Address of elements
  std::transform(vec_near_distances.begin(),vec_near_distances.end(),
                 back_inserter(vec_ref_near_distances),address<size_t>);
  std::transform(vec_far_distances.begin(),vec_far_distances.end(),
                 back_inserter(vec_ref_far_distances),address<size_t>);
  
  // Comparison functor
  dereferenced_less <size_t*> cmp_functor;

  // Sort intersections between spheres with different radii
  sort(vec_ref_near_distances.begin(), vec_ref_near_distances.end(),
       cmp_functor);
  sort(vec_ref_far_distances.begin(), vec_ref_far_distances.end(),
       cmp_functor);

#if _DEVELOPMENT


  std::vector<size_t*>::iterator it, next;
  size_t next_limit = SIZE_MAX;
  size_t last_limit = SIZE_MAX;
  

  // TEST this properly
  limits.push_back(data->m_minsample);
  bool bfirst = true;

  // Experimental
  if (!vec_ref_near_distances.empty()) {
    size_t bla = *vec_ref_near_distances.back();
    if (!(bla < data->m_minsample)) { // else not interesting

      it = find_if(vec_ref_near_distances.begin(),
		   vec_ref_near_distances.end(),
		   dereferenced_less_or_equal<size_t>(&data->m_minsample));

      if (it != vec_ref_near_distances.begin())
        it--;

      for (;it != vec_ref_near_distances.end();it++) {
        if (**it > data->m_minsample) {
          // Test for no repeated indices
          next = it;
          next++;
          if (next != vec_ref_near_distances.end())
            next_limit = **next;
          else
            next_limit = SIZE_MAX;

          if ((**it != last_limit) &&
              (**it != next_limit)) {
            if (bfirst) {
              order.push_back(it-vec_ref_near_distances.begin());
              bfirst = false;
            }
            order.push_back(it-vec_ref_near_distances.begin());
            limits.push_back(std::max(std::min(**it,
                                               max_index),data->m_minsample));
          }
        }
      }
    }
  }

  next_limit = SIZE_MAX;
  last_limit = SIZE_MAX;

  
  if (!vec_ref_far_distances.empty()) {
    if (!((*vec_ref_far_distances.back()) < data->m_minsample)) {
      // else not interesting

      std::vector<size_t*>::iterator it = 
				find_if(vec_ref_far_distances.begin(),
								vec_ref_far_distances.end(),
								dereferenced_less_or_equal<size_t>(&data->m_minsample));
			
      if (it != vec_ref_far_distances.begin())
        it--;
      for (;it != vec_ref_far_distances.end();it++) {
        if (**it > data->m_minsample) {
          // Test for no repeated indices
          next = it;
          next++;
          if (next != vec_ref_far_distances.end())
            next_limit = **next;
          else
            next_limit = SIZE_MAX;

          if ((**it != last_limit) &&
              (**it != next_limit)) {
            if (bfirst) {
              order.push_back(it-vec_ref_far_distances.begin());
              bfirst = false;
            }
            order.push_back(it-vec_ref_far_distances.begin());
            limits.push_back(std::max(std::min(**it,
                                               max_index),data->m_minsample));
          }
        }
      }      
    }
  }
#else
  // Works, but without minsample
  limits.push_back(0);

  bool bfirst = true;

  // Working reference
  if (!vec_ref_near_distances.empty()) {

    if ((*vec_ref_near_distances.back()) > 0) { // else not interesting

      // Works without minsample
      for (m=0;m<vec_ref_near_distances.size() ; m++) {
        if (*vec_ref_near_distances[m] > 0) {
          if (bfirst) {
            order.push_back(vec_ref_near_distances[m] - 
                            &vec_near_distances[0]);
            bfirst = false;
          }
          order.push_back(vec_ref_near_distances[m] -
                          &vec_near_distances[0]);
          // TEST: okay max(min(distance,max_index),min_index)
          limits.push_back(std::min(*vec_ref_near_distances[m],
				    max_index));
        }
      }
      
    }
  }
  
  if (!vec_ref_far_distances.empty()) {
    if ((*vec_ref_far_distances.back()) > 0) {
      for (m=0;m<vec_ref_far_distances.size() ; m++) {
        if (*vec_ref_far_distances[m] > 0) {
          if (bfirst) {
            order.push_back(vec_ref_far_distances[m] -  \
                            &vec_far_distances[0]);
            bfirst = false;
          }
          order.push_back(vec_ref_far_distances[m] -    \
                          &vec_far_distances[0]);
          limits.push_back(std::min(*vec_ref_far_distances[m],
				    max_index));
        }
      }
    }
  }
#endif

  // TODO: Fix this: Should be the one just further away
  if (bfirst) {
    // All points on the line lie further away than the furthest
    // distance, i.e. we use the apodization corresponding to being
    // further away than the largest distance
    if (!vec_ref_far_distances.empty()) {
      // Intersection found
      order.push_back(skip_indices+(vec_ref_far_distances.back() - \
                                    &vec_far_distances[0]));
    }
    else {
      // No intersection found
      order.push_back(max_dist_index);
    }
  }
  limits.push_back(max_index);

  /* Dynamic apodization is calculated using the distance from the
		 source or sink to the line rotated using the orientation given
		 using Euler angle
	euler_rot(apodization.data->m_apo_dir,this->data->m_direction,
						apodization.data->m_euler);
	*/

	/* Dynamic apodization is calculated using distance from source or
		 sink to the line passing though the focus point and oriented
		 according to the orientation given using Euler angles */
	//	basis_vector(apodization.data->m_apo_dir, apodization.data->m_euler,2);

	// TODO: Verify this
	euler[0] = apodization.data->m_aperture.data->m_euler[0+i_emission*3];
	euler[1] = apodization.data->m_aperture.data->m_euler[1+i_emission*3];
	euler[2] = apodization.data->m_aperture.data->m_euler[2+i_emission*3];
	basis_vector(apodization.data->m_apo_dir, euler,2);
}

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
