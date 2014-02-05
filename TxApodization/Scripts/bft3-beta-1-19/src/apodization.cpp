/*****************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Apodization class                                     *
 *                                                                            *
 * $Id: apodization.cpp,v 1.43 2011-07-17 21:49:10 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-07-17 21:49:10 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *   - Reference c from Line rather than aperture
 *   - Move validation from apodization_mex.cpp to Apodization
 *   -
 ****************************************************************************/

#include "aperture.h"
#include "apodization.h"

#ifdef _WIN32
# include <malloc.h>
#else
# include "mm_malloc.h"
#endif

#include <cstring>
#include <cmath>

// Used for delay calculation
#include <functional>
using std::greater;
using std::less;

//@{
/** Static variables */
size_t Apodization::DataRep::next_ID = 0;
//@}


Apodization::Apodization() : data(NULL) {}

Apodization::Apodization(const Apodization &other) {
  data = other.data;
  if (data)
    data->m_cnt++;
}

Apodization& Apodization::operator=(const Apodization &other) {
  if (other.data)
    other.data->m_cnt++;

  if (data)
    if (--data->m_cnt == 0) {
      delete data;
    }

  data = other.data;

  return *this;
}


Apodization::Apodization(const Aperture &aperture, const float_type* ref,
                         const float_type* distances, const size_t n_distances,
                         const float_type* values, const size_t n_elements)
  : data(NULL) {
  data = new Apodization::DataRep(ref, distances, n_distances,
                                  values, n_elements);
  data->m_aperture = aperture;
}

Apodization::~Apodization() {
  if (data)
    if (--data->m_cnt == 0)
      delete data;
}

void Apodization::clone() {
  data->m_cnt--;
  data = new DataRep(data->m_ref, data->m_distances, data->m_ndistances,
                     data->m_values, data->m_nelements);
}

Apodization Apodization::Clone() const {
  Apodization temp = Apodization(*this);

  // If data is used elsewhere
  if (data->m_cnt > 1) temp.clone();

  // Additional assigment of data (not covered by copy-constructor)
  temp.data->m_dynamic  = data->m_dynamic;
  temp.data->m_fixed    = data->m_fixed;
  temp.data->m_manual   = data->m_manual;
  temp.data->m_f        = data->m_f;
  temp.data->m_nactive_elements        = data->m_nactive_elements;
  temp.data->m_aperture = data->m_aperture;

  return temp;
}

// TODO: Remove this, see apodization_mex.cpp (Clone method)
Apodization* Apodization::Clone2() const {
  Apodization* temp = new Apodization(*this);

  // If data is used elsewhere
  if (data->m_cnt > 1) temp->clone();

  // Additional assigment of data (not covered by copy-constructor)
  temp->data->m_dynamic = data->m_dynamic;
  temp->data->m_fixed   = data->m_fixed;
  temp->data->m_manual  = data->m_manual;
  temp->data->m_f       = data->m_f;

  return temp;
}

bool Apodization::setAperture(const Aperture* aperture) {
  bool retval = false;
  if (aperture->data->m_npos == data->m_nelements) {
    data->m_aperture = *aperture;
    retval = true;
  }
  return retval;
}

Apodization::DataRep::DataRep() : m_ref(NULL),
                                  m_euler(NULL),
                                  m_distances(NULL),
                                  m_values(NULL),
                                  m_dynamic(false),
                                  m_manual(true),
                                  m_fixed(false),
                                  m_window((Apodization::Window_t)0),
                                  m_f(float_type(1.0)),
                                  m_id(next_ID++),
                                  m_cnt(1) {
  m_window_param = float_type(0.5);
}

Apodization::DataRep
::DataRep(const float_type* ref,
          const float_type* distances,
          const size_t n_distances,
          const float_type* values,
          const size_t n_elements) : m_ref(NULL),
                                     m_euler(NULL),
                                     m_distances(NULL),
                                     m_values(NULL),
                                     m_dynamic(false),
                                     m_manual(true),
                                     m_fixed(false),
                                     m_window((Apodization::Window_t)0),
                                     m_f(float_type(1.0)),
                                     m_id(next_ID++),
                                     m_cnt(1) {
  
  // TODO: Move validation from apodization_mex.cpp to Apodization
  m_ndistances = n_distances;
  m_nelements  = n_elements;
  m_nactive_elements = n_elements;
  m_dyn_size = float_type(0.05);
  m_window_param = float_type(0.5);

  m_ref       = (float_type*) _mm_malloc(3*sizeof(float_type),16);
  m_euler     = (float_type*) _mm_malloc(3*sizeof(float_type),16);
  m_distances = (float_type*) _mm_malloc(4*((n_distances+3)/4)*sizeof(float_type),16);
  m_values    = (float_type*) _mm_malloc(4*((n_elements*n_distances+3)/4)*sizeof(float_type),16);

  m_euler[0] = float_type(M_PI)/float_type(2.0);
  m_euler[1] = float_type(0.0); // rotation around y-axis if m_euler[0] = M_PI/2
  m_euler[2] = float_type(-M_PI)/float_type(2.0); // rotation in transducer plane (2D apodization)

  memcpy(m_ref,ref,3*sizeof(float_type));
  memcpy(m_distances,distances,m_ndistances*sizeof(float_type));
  memcpy(m_values,values,m_ndistances*m_nelements*sizeof(float_type));
}

// Destructor
Apodization::DataRep::~DataRep() {
  if (m_ref)
    _mm_free(m_ref);
  if (m_distances)
    _mm_free(m_distances);
  if (m_values)
    _mm_free(m_values);
  if (m_euler)
    _mm_free(m_euler);
}

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
