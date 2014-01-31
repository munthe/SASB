/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Aperture Class                                        *
 *                                                                            *
 * $Id: aperture.cpp,v 1.55 2011-08-04 18:29:16 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-08-04 18:29:16 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Inherit from base class, which defines xtors and assignments
 *  - Move initialization of focus to DataRep
 *  - Cannot clone aperture with more than 1 emission
 ****************************************************************************/

#include "aperture.h"

#include <cstring> // memcpy
#include <functional>
using std::greater;
using std::less;

#include "bfmath.h"

float_type Aperture::_fs = 100e6;
float_type Aperture::_f0 = 5e6;
float_type Aperture::_c  = 1540.0;

size_t Aperture::DataRep::next_ID = 0;

Aperture::Aperture() : data(NULL) {
}

Aperture::Aperture(const Aperture &other) {
  data = other.data;
  if (data)
    data->m_cnt++;
}

Aperture& Aperture::operator=(const Aperture &other) {
  if (other.data)
    other.data->m_cnt++;

  if (data)
    if (--data->m_cnt == 0) {
      delete data;
    }

  data = other.data;
  return *this;
}

Aperture::Aperture(const float_type* pos_vec, size_t length)  {
  data = new Aperture::DataRep(pos_vec, length, 1);
}

Aperture::Aperture(size_t n_elements, float_type pitch) {

  size_t i;

  float_type wx = float_type(n_elements-1)/2;
           
  float_type* pos_vec =
    (float_type*) _mm_malloc(4*((n_elements*3+3)/4)*sizeof(float_type),16);
  
  for (i=0;i<n_elements;i++) {
    pos_vec[i] = (float_type(i) - wx)*pitch;
  }
  for (i=n_elements;i<3*n_elements;i++) {
    pos_vec[i] = float_type(0.0);
  }
  
  data = new Aperture::DataRep(pos_vec,3*n_elements,1,Aperture::linear_array);
  
  data->m_pitch = pitch;
	data->m_radius = float_type(0);

  _mm_free(pos_vec);
}

Aperture::Aperture(size_t n_elements, float_type width, float_type kerf,
									 float_type radius) {

  size_t i;

  float_type* pos_vec =
    (float_type*) _mm_malloc(4*((n_elements*3+3)/4)*sizeof(float_type),16);

  float theta_pitch = asin (width/radius) + asin (kerf/radius);

  float theta;

  float_type wtheta = float_type(n_elements-1)/2;

  for (i=0;i<n_elements;i++) {
		theta = (float_type(i) - wtheta)*(theta_pitch);
    pos_vec[i]              = radius * sin(theta);
		pos_vec[i+n_elements]   = 0;
		pos_vec[i+2*n_elements] = -radius * (float_type(1)-cos(theta));
	}

  data = new Aperture::DataRep(pos_vec,3*n_elements,1,Aperture::convex_array);
  
  data->m_pitch = 2*radius*sin(theta_pitch/2);
	data->m_radius = radius;
  _mm_free(pos_vec);
}

Aperture::Aperture(size_t n_xelements, size_t n_yelements,
									 float_type xpitch, float_type ypitch) {

  size_t i,j;

  float_type wx = float_type(n_xelements-1)/2;
  float_type wy = float_type(n_yelements-1)/2;
           
  float_type* pos_vec =
    (float_type*) _mm_malloc(n_xelements*n_yelements*3*sizeof(float_type),16);
  
  for (j=0;j<n_yelements;j++) {
		for (i=0;i<n_xelements;i++) {
			pos_vec[j*n_xelements + i]     = (float_type(i) - wx)*xpitch;
			pos_vec[j*n_xelements + i + 1] = (float_type(j) - wy)*ypitch;
			pos_vec[j*n_xelements + i + 2] = float_type(0.0);
		}
  }
  
  data = new Aperture::DataRep(pos_vec,3*n_xelements*n_yelements,1,
															 Aperture::two_dim_array);
  
	/* TODO: Support 2D pitch for fancy apodizations */
  data->m_pitch = xpitch;
	data->m_radius = float_type(0);

  _mm_free(pos_vec);
}

Aperture::~Aperture() {
  if (data)
    if (--data->m_cnt == 0)
      delete data;
}

// TODO: Move to DataRep and inline in header
bool Aperture::setPositions(const float_type* pos, const size_t npos) {
  bool retval = false;
  if ((data->m_npos*3) == npos) {
    memcpy(&data->m_pos->m_data[0][0],pos,npos*sizeof(float_type));
    retval = true;
  }
  return retval;
}

void Aperture::clone() {
  data->m_cnt--;
  data = new DataRep(&data->m_pos->m_data[0][0],data->m_npos*3, data->m_nemissions);
}

Aperture* Aperture::Clone() const {
  size_t i;
  Aperture *temp = new Aperture(*this);

  // If data is used elsewhere
  if (data->m_cnt > 1)
    temp->clone();

  // Check operator=

  // Additional assigment of data (not covered by copy-constructor)
  if (data->m_focus) {
    temp->data->m_focus = (float_type*) _mm_malloc(sizeof(float_type)*4,16);
    memcpy(temp->data->m_focus, data->m_focus, sizeof(float)*4);
  }
  if (data->m_delays) {
    temp->data->m_delays = 
      (float_type*) _mm_malloc(sizeof(float_type)*data->m_npos,16);
    memcpy(temp->data->m_delays,data->m_delays,
           sizeof(float)*data->m_npos);
  }

	/* TODO: Use memcpy */
  for (i=0; i < 3*data->m_nemissions;i++) {
    temp->data->m_center_focus[0+i*3] = data->m_center_focus[0+i*3];
    temp->data->m_center_focus[1+i*3] = data->m_center_focus[1+i*3];
    temp->data->m_center_focus[2+i*3] = data->m_center_focus[2+i*3];
    temp->data->m_euler[0+i*3] = data->m_euler[0+i*3];
    temp->data->m_euler[1+i*3] = data->m_euler[1+i*3];
    temp->data->m_euler[2+i*3] = data->m_euler[2+i*3];
  }
  temp->data->m_pitch  = data->m_pitch;
  temp->data->m_radius = data->m_radius;
  temp->data->m_type   = data->m_type;

  return temp;
}

bool Aperture::getFocusDelays(float_type* focus_delays) {
  bool retval = false;
  size_t i;
  if (this->data->m_focus) {
    for (i=0;i<this->data->m_npos;i++) {
      focus_delays[i] = sqrt(SQUARE((*this)(0,i)-this->data->m_focus[0]) +
                             SQUARE((*this)(1,i)-this->data->m_focus[1]) +
                             SQUARE((*this)(2,i)-this->data->m_focus[2]));
      focus_delays[i] = focus_delays[i] / (*this->data->m_c);
    }

    float_type center_delay =
      *std::max_element(focus_delays, 
                        focus_delays + this->data->m_npos,less<float_type>());
    for (i=0;i<this->data->m_npos;i++) {
      focus_delays[i] = center_delay - focus_delays[i];
    }
    retval = true;
  }
  return retval;
}

// TODO: Use STL containers instead
bool Aperture::setCenterFocus(const float_type* center_focus,
                              size_t n_emissions) {
  // Use content from aperture_mex.cpp
  size_t i;
  if (data->m_nemissions != n_emissions) {
    _mm_free(data->m_center_focus);
    data->m_center_focus =
      (float_type*) _mm_malloc(3*n_emissions*sizeof(float_type),16);

		// Old Euler angles
		float_type* p_euler = data->m_euler;

		// Re-allocate Euler angles and copy old ones
		if (data->m_euler) {
			data->m_euler = (float_type*) _mm_malloc(3*n_emissions*sizeof(float_type),16);
			memset(data->m_euler,0,3*n_emissions*sizeof(float_type));

			for (i=0;i<std::min(data->m_nemissions,n_emissions);i++) {
				data->m_euler[0+i]            = p_euler[0+i];
				data->m_euler[n_emissions+i]  = p_euler[data->m_nemissions+i];
				data->m_euler[2*n_emissions+i]= p_euler[2*data->m_nemissions+i];
			}
			_mm_free(p_euler);
		}
    data->m_nemissions = n_emissions;
  }
  for (i=0 ; i < data->m_nemissions;i++) {
    data->m_center_focus[0+i*3] = center_focus[i];
    data->m_center_focus[1+i*3] =
      center_focus[data->m_nemissions+i];
    data->m_center_focus[2+i*3] =
      center_focus[2*data->m_nemissions+i];
  }
  
  return true;
}

bool Aperture::setOrientation(const float_type* euler,
                              size_t n_emissions) {
  // Use content from aperture_mex.cpp
  size_t i;
  if (data->m_nemissions != n_emissions) {
    _mm_free(data->m_euler);
    data->m_euler =
      (float_type*) _mm_malloc(3*n_emissions*sizeof(float_type),16);
		float_type* p_center_focus = data->m_center_focus;

		// Re-allocate center_focus
		if (data->m_center_focus) {
			data->m_center_focus = (float_type*) _mm_malloc(3*n_emissions*sizeof(float_type),16);
			memset(data->m_center_focus,0,3*n_emissions*sizeof(float_type));

			for (i=0;i<std::min(data->m_nemissions,n_emissions);i++) {
				data->m_center_focus[0+i]            = p_center_focus[0+i];
				data->m_center_focus[n_emissions+i]  = p_center_focus[data->m_nemissions+i];
				data->m_center_focus[2*n_emissions+i]= p_center_focus[2*data->m_nemissions+i];
			}
			_mm_free(p_center_focus);
		}

    data->m_nemissions = n_emissions;
  }

  for (i=0 ; i < data->m_nemissions;i++) {
    data->m_euler[0+i*3] = euler[i];
    data->m_euler[1+i*3] =
      euler[data->m_nemissions+i];
    data->m_euler[2+i*3] =
      euler[2*data->m_nemissions+i];
  }
  
  return true;
}

// Verify, want to copy if referenced elsewhere (may be not)
bool Aperture::setFocus(const float_type* focus) {

  if (focus) {
    if (data->m_focus)
      memcpy(data->m_focus, focus, 3*sizeof(float_type));
    else {
      data->m_focus = (float_type*) _mm_malloc(sizeof(float_type)*4,16);
      memcpy(data->m_focus, focus, 3*sizeof(float_type));    
    }
  }
  else {
    if (data->m_focus) {
      _mm_free(data->m_focus);
      data->m_focus = NULL;
    }
  }
  return true;
}

bool Aperture::setDelays(const float_type* delays) {
  size_t i;

  if (delays) {
    if (!data->m_delays) {
      data->m_delays =
        (float_type*) _mm_malloc(sizeof(float_type)*n_elements(),
                                 16);
    }
    for (i=0 ; i < n_elements() ; i++)
      data->m_delays[i] = delays[i];
  }
  else {
    if (data->m_delays) {
      _mm_free(data->m_delays);
      data->m_delays = NULL;
    }
  }
  return true;
}

Aperture::DataRep::DataRep() : m_pos(NULL), m_center_focus(NULL), m_euler(NULL), m_focus(NULL),
                               m_delays(NULL), m_fs(NULL), m_f0(NULL),
                               m_c(NULL), m_id(next_ID++), m_type(Aperture::custom), m_ppwave(false), m_cnt(1)
{
  m_fs = &Aperture::_fs;
  m_f0 = &Aperture::_f0;
  m_c  = &Aperture::_c;

	m_nemissions = 1;
  m_center_focus = (float_type*) _mm_malloc(m_nemissions*3*sizeof(float_type),16);
  m_euler        = (float_type*) _mm_malloc(m_nemissions*3*sizeof(float_type),16);
}
  
Aperture::DataRep::DataRep(const float_type* pos_vec,
                           size_t length, size_t n_emissions, Aperture::Type type) : 
  m_pos(NULL), m_center_focus(NULL), m_euler(NULL), m_focus(NULL), m_delays(NULL),
  m_fs(NULL), m_f0(NULL), m_c(NULL), m_id(next_ID++), m_type(type), m_ppwave(false), m_cnt(1)
{
  
  m_fs = &Aperture::_fs;
  m_f0 = &Aperture::_f0;
  m_c  = &Aperture::_c;
  
  m_npos = length/3; 
  m_nemissions = n_emissions;

  m_pos = new Aperture::PMatrix<float_type>(m_npos,3);

  m_center_focus = (float_type*) _mm_malloc(n_emissions*3*sizeof(float_type),16);
  m_euler        = (float_type*) _mm_malloc(n_emissions*3*sizeof(float_type),16);

	memset(m_center_focus,0,n_emissions*3*sizeof(float_type));
	memset(m_euler,0,n_emissions*3*sizeof(float_type));

  memcpy(&m_pos->m_data[0][0],pos_vec,length*sizeof(float_type));

  if (m_type==Aperture::custom) {
    if (m_npos > 1)  {
      this->m_center_focus[0] = (pos_vec[m_npos-1]+pos_vec[0])/2;
      this->m_center_focus[1] = (pos_vec[2*m_npos-1]+pos_vec[m_npos])/2;
      this->m_center_focus[2] = (pos_vec[3*m_npos-1]+pos_vec[2*m_npos])/2;

      this->m_pitch = sqrt(SQUARE(pos_vec[0]-pos_vec[1]) + 
                           SQUARE(pos_vec[m_npos]-pos_vec[m_npos+1]) + 
                           SQUARE(pos_vec[2*m_npos]-pos_vec[2*m_npos+1]));
    }
    else {
      this->m_center_focus[0] = pos_vec[0];
      this->m_center_focus[1] = pos_vec[1];
      this->m_center_focus[2] = pos_vec[2];
      this->m_pitch = float_type(0.0);
    }
  }
  else {
    this->m_center_focus[0] = (pos_vec[m_npos-1]+pos_vec[0])/2;
    this->m_center_focus[1] = float_type(0.0);
    this->m_center_focus[2] = float_type(0.0);
  }
}

Aperture::DataRep::~DataRep() {
  if (m_pos)
    delete m_pos;
  if (m_focus)
    _mm_free(m_focus);
  if (m_center_focus) {
    _mm_free(m_center_focus);
  }
  if (m_euler) {
    _mm_free(m_euler);
  }
  if (m_delays)
    _mm_free(m_delays);
}

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
