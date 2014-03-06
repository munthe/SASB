/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Line Class                                            *
 *                                                                            *
 * $Id: line.h,v 1.42 2011-07-19 20:17:24 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-19 20:17:24 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 * $Log: line.h,v $
 * Revision 1.42  2011-07-19 20:17:24  jmh
 * Moved orientation to aperture class
 *
 * Revision 1.41  2011/04/16 13:05:26  jmh
 * *** empty log message ***
 *
 * Revision 1.40  2011/04/06 22:12:43  jmh
 * *** empty log message ***
 *
 * Revision 1.39  2011/04/06 13:52:04  jmh
 * Error when using manual apodization
 *
 * Revision 1.36  2010/09/02 20:15:28  jmh
 * Removed line delays
 *
 * Revision 1.35  2010/09/02 19:42:58  jmh
 * Line delays
 *
 * Revision 1.33  2010/07/17 15:00:03  jmh
 * line apodization property
 *
 * Revision 1.29  2010/04/09 14:15:39  jmh
 * Problem with fixed emit apodization
 *
 * Revision 1.24  2010/03/30 19:34:15  jmh
 * VS XMT and RCV - slow but it works
 *
 * Revision 1.22  2010/03/27 19:54:12  jmh
 * Main thread initializes interpolators
 *
 * Revision 1.21  2010/03/20 13:52:12  jmh
 * 2nd attempt dyn emit
 *
 * Revision 1.19  2010/03/19 19:57:07  jmh
 * First attempt dyn xmt apo
 *
 * Revision 1.18  2010/03/19 17:52:26  jmh
 * AV when beamforming t1.m
 *
 * Revision 1.17  2010/03/14 03:20:08  jmh
 * TODO: BLine(Apodization,....)
 *
 * Revision 1.6  2009/11/04 13:23:21  jmh
 * Removed corrupted double linked list, now just wrong answers
 *
 * Revision 1.5  2009/11/04 09:45:01  jmh
 * Corrupted double linked list
 *
 * Revision 1.4  2009/10/27 14:33:08  jmh
 * single precision slow attempt
 *
 * Revision 1.2  2009/10/12 21:05:53  jmh
 * Fixed DataRep
 *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 ****************************************************************************/

#ifndef LINE_H
#define LINE_H

#include "common.h"

#include <vector>
using std::vector;

class Apodization;

template <typename T>
T* address(T& t) { return &t;}

class Line {
 public:

  template <class T> class dereferenced_less {
  public:
    bool operator() (const T lhs, const T rhs) const {return *lhs < *rhs;}
  };

  template<class T> class dereferenced_less_or_equal {
  public:
    dereferenced_less_or_equal(const T* p) : p_(p) { }
    bool operator()(const T* p2) const {return (*p_ < *p2); }
  private:
    const T* p_;
  };
  
  class sample_apodization_create {
  public:
  sample_apodization_create(int apo) : m_apo(apo){};
    std::pair<size_t, int> operator()(size_t t) {
      return std::make_pair(t,m_apo);
    }
  private:
    int m_apo;
  };
  
/*
  class find_apo : public std::unary_function<std::pair<size_t, int>, bool> {
  public:
  find_apo(const int apo) : m_apo(apo) {}
    
    bool inline operator()(std::pair<size_t, int>& sample_apo_pair) const {
      return (sample_apo_pair.second == m_apo);
    }
  private:
    int m_apo;
  };
*/

  /** 
   * Default constructor
   * 
   */
  Line();

  /** 
   * Constructor:
   * 
   * @param origin 
   * @param direction 
   * @param dr 
   * @param length 
   */
  Line(const float_type* origin, const float_type* direction,
       const float_type dr, const float_type length);

  /** 
   * Constructor: Sets up reference for apodization
   * 
   * @param xmt_apodization 
   * @param rcv_apodization 
   * @param origin 
   * @param direction 
   * @param dr 
   * @param length 
   */
  Line(const Apodization& xmt_apodization,
       const Apodization& rcv_apodization,
       const float_type* origin,
       const float_type* direction, const float_type dr,
       const float_type length);

  /** 
   * Copy-constructor
   * 
   * @param other 
   */
  Line(const Line &other);

  /** 
   * 
   * 
   * @param other 
   * 
   * @return 
   */
  Line& operator=(const Line &other);
  /** 
   * Destructor:
   * 
   */
  ~Line();

  void findMinSample(const float_type delay, const Aperture& aperture);

  /** 
   * Working reference for setting up pixel indices for a given apodization
   * 
   * @param aperture 
   * @param apodization 
   * @param delay 
   */
  void setAnyApodization(const Apodization& apodization,
			 vector<size_t>& order,
			 vector<size_t>& limits,
			 size_t i_emission=0);
  /** 
   * Find the samples crossing the virtual plane
   * 
   * @param apodization 
   * @param vp_limits 
   * @param order 
   */
  void findVirtualPlane(const Apodization& apodization,
			vector<size_t>& vp_limits,
			vector<size_t>& order);

  /** 
   * Initialize apodization
   * 
   * @param xmt_apodization 
   * @param rcv_apodization 
   * @param trans_no 
   */
  void setApodization2(const Apodization& xmt_apodization,
		       const Apodization& rcv_apodization,
		       const size_t trans_no);

  /** 
   * Compute apodization value for either transmit or receive. It is
   * assumed that the setApodization function has been called
   * previously
   * 
   * @param apo_values 
   * @param direction 
   * @param trans_no 
   * 
   * @return 
   */
  bool apodize(float_type* apo_values, int direction, size_t trans_no=0);


  /* TODO: Replace with references & */
  inline float_type length() const { return data->m_length; }
  inline float_type dr() const { return data->m_dr; }

  inline float_type dx() const { return data->m_dr*data->m_direction[0]; }
  inline float_type dy() const { return data->m_dr*data->m_direction[1]; }
  inline float_type dz() const { return data->m_dr*data->m_direction[2]; }

  inline float_type operator()(size_t dim, size_t dist) {
    return this->origin(dim) + data->m_dr*data->m_direction[dim]*dist;
  }

  inline float_type direction(size_t dim) const { return data->m_direction[dim]; }
  inline float_type origin(size_t dim) const { return data->m_origin[dim]; }

  class DataRep {
  public:   
    /** 
     * 
     * 
     * 
     * @return 
     */
    DataRep();
    /** 
     * 
     * 
     * @param origin 
     * @param direction 
     * @param dr 
     * @param length 
     * 
     * @return 
     */
    DataRep(const float_type* origin, const float_type* direction,
	    const float_type dr, const float_type length);

    /** 
     * 
     * 
     * @param apodization 
     * @param origin 
     * @param direction 
     * @param dr 
     * @param length 
     * 
     * @return 
     */
    DataRep(const Apodization& xmt_apodization,
	    const Apodization& rcv_apodization,
	    const float_type* origin, const float_type* direction,
	    const float_type dr, const float_type length);
    
    /** 
     * 
     * 
     * 
     * @return 
     */
    ~DataRep();

    // Usage count
  public:

    size_t m_cnt;

    float_type* m_origin;
    float_type* m_direction;
    float_type m_dr;
    float_type m_length;

    size_t m_npos;
    size_t m_minsample;

    vector<size_t> m_rcv_apo_order; // Index for apodization
    vector<size_t> m_xmt_apo_order; // Index for apodization
    vector<size_t> m_rcv_vs_order;  // Index for VS calculation
    vector<size_t> m_xmt_vs_order;  // Index for VS calculation
    vector<size_t> m_limits; // Index for apodization

    Apodization m_xmt_apodization;
    Apodization m_rcv_apodization;

  };

  DataRep* data;
};

#endif
