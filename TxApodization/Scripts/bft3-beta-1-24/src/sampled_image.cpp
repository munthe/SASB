/******************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Sampled Image Class                                   *
 *                                                                            *
 * $Id: sampled_image.cpp,v 1.40 2011-08-05 13:23:43 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-08-05 13:23:43 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 * - Verify cloning of apodization
 * - Allow for manual apodization, Aperture centered around BF line
 *   (synthesized aperture), fixed by specifying distance as well as a
 *   number of windows (orthogonal to BF line)
 * - How to specify which apodization goes to which points?
 * - 3D points 
 * - Overload new and delete or use explicit
 *   ctor / dtor of SampleInterpolate1D
 * - Change to support both multiple angles and emissions (now only angles)
 ****************************************************************************/
#include "aperture.h"
#include "apodization.h"
#include "sampled_image.h"

#ifdef __GNUC__
# ifndef HAVE_THREAD
#  define HAVE_THREAD
# endif
#endif

// Debugging
#include <cstdio>
#include <cstdlib>

// Error handling
#include <cerrno>
#include <cstring> // strerror

#include "sample_interpolate.h"

#include "mex.h"

#ifdef __linux__
static pid_t gettid( void ) {
  pid_t pid;
  CallErr(pid = syscall, (__NR_gettid));
  return pid;
}
#else
#endif

typedef struct thread_arg {
  size_t start_x;
  size_t end_x;
  float_type* image;
  const float_type* rf_data;
  const float_type* delay;
  Apodization xmt_apo;
  Apodization rcv_apo;
  const float_type* angles;
  size_t n_angles;
  float_type* xs;
  float_type* zs;
  float_type* ys;
  size_t n_rf_samples;
  const uint32_t* trans_no;
  size_t thread_id;
  int cpu_id;
} thread_arg_t;

// Hack, remove after IEEE
inline void convex_apo_ref(const float_type* pixel, const float_type angle,
													 const float_type radius, float_type* ref) __attribute__ ((always_inline));

inline void convex_apo_ref(const float_type* pixel, const float_type angle,
													 const float_type radius, float_type* ref) {
	float_type p1[2];
	float_type p2[2];
	float_type dx, dy, dr, D, det;
	float_type tmp1 = 0;
	float_type tmp2 = 0;

	float_type _angle = angle - atan2(pixel[0],pixel[2]+radius);

	p1[0] = pixel[0];
	p1[1] = pixel[2];
	p2[0] = p1[0] - sin(_angle); // Verify sign
	p2[1] = p1[1] + cos(_angle);
	
	// Subtract origin vector
	p1[0] -= 0;
	p2[0] -= 0;
	p1[1] -= -radius;
	p2[1] -= -radius;

	dx = p2[0] - p1[0];
	dy = p2[1] - p1[1];
	dr = sqrt(SQUARE(dx)+SQUARE(dy));
	D = p1[0]*p2[1] - p2[0]*p1[1];
	det = SQUARE(radius)*SQUARE(dr)-SQUARE(D);

	if (det < 0) {
		ref[0] = float_type(0);
		ref[1] = float_type(0);
		ref[2] = float_type(0);
	}
	else {
		tmp1 = -D*dx+fabs(dy)*sqrt(det)/(SQUARE(dr)) + (-radius); // Added origin vector
		tmp2 = -D*dx-fabs(dy)*sqrt(det)/(SQUARE(dr)) + (-radius); // Added origin vector
		if (fabs(tmp1) > fabs(tmp2)) {
			ref[2] = tmp2;
			ref[1] = float_type(0);
			ref[0] = D*dy-copysign(1.0,dy)*dx*sqrt(det)/(SQUARE(dr));
		}
		else {
			ref[2] = tmp1;
			ref[1] = float_type(0);
			ref[0] = D*dy+copysign(1.0,dy)*dx*sqrt(det)/(SQUARE(dr));
		}
	}
}

// Template instantiation
template class pthread_launcher<SampledImage,
                                &SampledImage::beamform_threaded>;

size_t SampledImage::nthreads = 1;

SampledImage::SampledImage() {
}

SampledImage::SampledImage(const Aperture& Th_t, const Aperture& Th_r,
													 const Apodization& Ah_t, const Apodization& Ah_r,
													 sampled_image_t im)
{

  interp_type = SampledImage::linear; // Default interpolation method
  m_xmt_aperture = Aperture(Th_t);
  m_rcv_aperture = Aperture(Th_r);
  m_xmt_apodization = Apodization(Ah_t);
  m_rcv_apodization = Apodization(Ah_r);

  im_geom = im;
}

bool SampledImage::Beamform(float_type* image, const float_type* rf_data,
                            const float_type* delay, const size_t n_rf_samples,
                            const uint32_t* trans_no,
                            size_t n_angles, const float_type* angles) {

  size_t i;
  int nproc;

#ifdef __linux__
  cpu_set_t cpuset;
#endif

#if defined(HAVE_PTHREAD_H)
  pthread_attr_t attr;
  pthread_t threads[N_MAX_THREADS];
#else
  unsigned int threadID;
  HANDLE threads[N_MAX_THREADS];
#endif

  thread_arg_t threadarg[N_MAX_THREADS];

  pthread_launcher<SampledImage,
    &SampledImage::beamform_threaded> launcher[N_MAX_THREADS];

#if defined(__linux__)
  CPU_ZERO( &cpuset );
  CallErr(sched_getaffinity,
          (gettid(), sizeof( cpu_set_t ), &cpuset));
  nproc = CPU_COUNT(&cpuset);
#elif defined(_MSC_VER)
  SYSTEM_INFO info;
  GetSystemInfo(&info);
  nproc = info.dwNumberOfProcessors;
#endif

#ifdef HAVE_PTHREAD_H
  CallErr(pthread_attr_init,(&attr));
  CallErr(pthread_attr_setdetachstate, (&attr, PTHREAD_CREATE_JOINABLE));
#endif

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

  size_t i_x, i_z, i_y;

  // Coordinate position (TODO: Calculate on-the-fly)
  float_type wx =
    (float_type(im_geom.nx)-float_type(1.0))/float_type(2.0) + im_geom.offset_x;
  float_type wz =
    (float_type(im_geom.nz)-float_type(1.0))/float_type(2.0) + im_geom.offset_z;
  float_type wy =
    (float_type(im_geom.ny)-float_type(1.0))/float_type(2.0) + im_geom.offset_y;

  float_type* xs =
    (float_type*) _mm_malloc( (size_t) ((im_geom.nx) * sizeof(float_type)),16);
  float_type* zs =
    (float_type*) _mm_malloc( (size_t) ((im_geom.nz) * sizeof(float_type)),16);  
  float_type* ys =
    (float_type*) _mm_malloc( (size_t) ((im_geom.ny) * sizeof(float_type)),16);

  for (i_x=0;i_x<im_geom.nx;i_x++)
    xs[i_x] = im_geom.dx * (((float_type) i_x)-wx);
  for (i_z=0;i_z<im_geom.nz;i_z++)
    zs[i_z] = im_geom.dz * (((float_type) i_z)-wz);
  for (i_y=0;i_y<im_geom.ny;i_y++)
    ys[i_y] = im_geom.dy * (((float_type) i_y)-wy);

  // Clone apodization objects
  Apodization* xmt_apodization =
    (Apodization*)_mm_malloc(nthreads*sizeof(Apodization),16);
  Apodization* rcv_apodization =
    (Apodization*)_mm_malloc(nthreads*sizeof(Apodization),16);

  for (i=0;i<nthreads;i++) {
    new (xmt_apodization+i) Apodization();
    new (rcv_apodization+i) Apodization();
    xmt_apodization[i] = m_xmt_apodization.Clone();
    rcv_apodization[i] = m_rcv_apodization.Clone();
  }

  for (i=0;i<nthreads;i++) {
    // Populate structs
    threadarg[i].start_x = 0+i*(im_geom.nx/nthreads);
    threadarg[i].end_x   = (im_geom.nx/nthreads)+i*(im_geom.nx/nthreads);
    threadarg[i].image = image;
    threadarg[i].xs = xs;
    threadarg[i].zs = zs;
    threadarg[i].ys = ys;
    threadarg[i].rf_data = rf_data;
    threadarg[i].delay = delay;
    threadarg[i].xmt_apo = xmt_apodization[i];
    threadarg[i].rcv_apo = rcv_apodization[i];
    threadarg[i].n_angles = n_angles;
    threadarg[i].angles = angles;
    threadarg[i].n_rf_samples = n_rf_samples;
    threadarg[i].trans_no = trans_no;
    threadarg[i].thread_id = i;
    threadarg[i].cpu_id = ((int) i) % nproc;
    if (i==(nthreads-1)) {
      threadarg[i].end_x = im_geom.nx;
    }
    launcher[i] =
      pthread_launcher<SampledImage,
      &SampledImage::beamform_threaded>(this,
																				&threadarg[i]);
  }

  // Reset output image
  memset(image,0,im_geom.nx*im_geom.nz*sizeof(float_type));

  for (i=0;i<nthreads;i++) {
#if defined(HAVE_PTHREAD_H)
    CallErr(pthread_create,
            (&threads[i], &attr,
             launch_member_function<pthread_launcher<SampledImage,
             &SampledImage::beamform_threaded> >, 
             &launcher[i]));
#elif defined(_MSC_VER)
    threads[i] =
      (HANDLE)_beginthreadex(NULL, 0,
                             launch_member_function<pthread_launcher<SampledImage,
                             &SampledImage::beamform_threaded> >,
                             &launcher[i], 0, &threadID );
#endif
  }
  for (i = 0; i < nthreads; i++) {
#if defined(HAVE_PTHREAD_H)
    CallErr(pthread_join,(threads[i],NULL));
#elif defined(_MSC_VER)
    WaitForSingleObject(threads[i], INFINITE );
#endif
  }
  
  for (i=0;i<nthreads;i++) {
    xmt_apodization[i].~Apodization();
    rcv_apodization[i].~Apodization();
  }
  _mm_free(xmt_apodization);
  _mm_free(rcv_apodization);

  _mm_free(xs);
  _mm_free(zs);
  _mm_free(ys);
    
#if defined(HAVE_PTHREAD_H)
  CallErr(pthread_attr_destroy,(&attr));
#endif
  
  return true;
}

// TODO: Speciale according to apodization
#if defined(HAVE_PTHREAD_H)
void* SampledImage::beamform_threaded(void *ptarg)
#else
  unsigned __stdcall SampledImage::beamform_threaded(void *ptarg)
#endif
{
  thread_arg_t* arg = (thread_arg_t*) ptarg;
  
  size_t i_xmt, i_x, i_z, i_rcv, i_angle;
  size_t i_y;

  size_t i;
  const float_type* rf_data      = arg->rf_data;

  // Local variables, TODO: Make some const
  const size_t x_start           = arg->start_x;
  const size_t x_end             = arg->end_x;

  float_type* image              = arg->image;
  const size_t nz                = im_geom.nz;

  Apodization xmt_apo = arg->xmt_apo;
  Apodization rcv_apo = arg->rcv_apo;

  float_type angle;
  const float_type* angles     = arg->angles;
  size_t n_angles              = arg->n_angles;
  const float_type* delay    = arg->delay;
  const size_t n_rf_samples  = arg->n_rf_samples;

  // TODO: Compute on-the-fly
  const float_type* xs       = arg->xs;
  const float_type* zs       = arg->zs;
  const float_type* ys       = arg->ys;

  float_type dist_xmt, dist_rcv, result, result1;
  result  = float_type(0);
	result1 = float_type(0);

  float_type findex;
  size_t trans_no;

  // Speed of sound and sample rate
  const float_type c  = *m_xmt_aperture.data->m_c;
  const float_type fs = *m_xmt_aperture.data->m_fs;

  const size_t n_rcv_elements = m_rcv_aperture.n_elements();

  const float_type rcv_ap_size = m_rcv_aperture.data->m_pitch*float_type(rcv_apo.data->m_nactive_elements);

	//	mexPrintf("rcv_ap_size: %f\n",rcv_ap_size);

	// IEEE Hack
	m_xmt_aperture.data->m_radius = m_rcv_aperture.data->m_radius;
	m_xmt_aperture.data->m_pitch  = m_rcv_aperture.data->m_pitch;

	const float_type xmt_ap_size = m_xmt_aperture.data->m_pitch*float_type(xmt_apo.data->m_nactive_elements);

  //  mexPrintf("n_active_elements: %d\n",xmt_apo.data->m_nactive_elements);
  ALIGN16_BEGIN float_type xyz_t[3] ALIGN16_END;

  ALIGN16_BEGIN float_type pixel[3] ALIGN16_END;


  for (i_xmt=0; i_xmt < m_xmt_aperture.data->m_nemissions ; i_xmt++) {

    trans_no = (arg->trans_no[i_xmt]-1);

		/* TODO: Share interpolator between threads, tricky */
    SampleInterpolate<size_t,float_type>** interpl = NULL;

    interpl =
      new SampleInterpolate<size_t,float_type>*[n_rcv_elements];
    for (i=0; i < n_rcv_elements ; i++) {
      interpl[i] =
        new SampleInterpolate<size_t,
															float_type>(&rf_data[i*n_rf_samples/*+ i_xmt*n_rf_samples*n_rcv_elements*/], n_rf_samples);
      interpl[i]->setMethod(interp_type);
    }
    xyz_t[0] = m_xmt_aperture.data->m_center_focus[0+i_xmt*3];
    xyz_t[1] = m_xmt_aperture.data->m_center_focus[1+i_xmt*3];
    xyz_t[2] = m_xmt_aperture.data->m_center_focus[2+i_xmt*3];

    for (i_x=x_start;i_x<x_end;i_x++) {
      pixel[0] = xs[i_x];
      for (i_y=0;i_y<1;i_y++) {
        pixel[1] = ys[i_y]; // TODO: Support y-coordinates
        for (i_z=0;i_z<nz;i_z++) {

          // TODO: Move apodization innermost
          pixel[2] = zs[i_z];

          // XMT distance
          dist_xmt = sqrt(SQUARE(xs[i_x] - xyz_t[0]) + 
                          SQUARE(zs[i_z] - xyz_t[2]) +
                          SQUARE(ys[i_y] - xyz_t[1]));
          
          dist_xmt = fs*dist_xmt/c;

          for (i_rcv=0 ; i_rcv< n_rcv_elements ; i_rcv++) {
            dist_rcv = sqrt(SQUARE(xs[i_x] - m_rcv_aperture(0,i_rcv)) +
                            SQUARE(zs[i_z] - m_rcv_aperture(2,i_rcv)) +
                            SQUARE(ys[i_y] - m_rcv_aperture(1,i_rcv)));

            dist_rcv = fs*dist_rcv/c;
            
            findex = dist_rcv + dist_xmt;
            
            findex = findex - delay[i_xmt]*fs;

            if ((findex > float_type(0.0)) && (findex < (n_rf_samples-1))) {
              
              result1 = (*interpl[i_rcv])(findex);

              for (i_angle=0;i_angle<n_angles;i_angle++) {
								result = result1;
                angle = angles[i_angle];

                // TODO: Support two different apertures, this assumes
                // aperture positioned in xy-plane
                float_type recv_offset = (zs[i_z]*tan(angle) + xs[i_x]);
                
                // Receive apodization
								if (m_rcv_aperture.data->m_type==Aperture::linear_array) {
									rcv_apo.data->m_ref[0] = recv_offset;
									rcv_apo.data->m_ref[1] = 0.0;
									rcv_apo.data->m_ref[2] = 0.0;
								}
								else if (m_rcv_aperture.data->m_type==Aperture::convex_array) {
									convex_apo_ref(pixel, angle, m_rcv_aperture.data->m_radius,
																 rcv_apo.data->m_ref);
								}
								else {
									rcv_apo.data->m_ref[0] = recv_offset;
									rcv_apo.data->m_ref[1] = 0.0;
									rcv_apo.data->m_ref[2] = 0.0;
								}
								/* Added because of introduction of m_apo_dir for
									 transmit apodization, TODO: Verify */
								rcv_apo.data->m_apo_dir[0] = -sin(angle);
								rcv_apo.data->m_apo_dir[1] = 0;
								rcv_apo.data->m_apo_dir[2] = cos(angle);
								
								if (m_xmt_aperture.data->m_type==Aperture::linear_array) {
									xmt_apo.data->m_ref[0] = recv_offset;
									xmt_apo.data->m_ref[1] = 0.0;
									xmt_apo.data->m_ref[2] = 0.0;
								}
								else if (m_xmt_aperture.data->m_type==Aperture::convex_array) {
									convex_apo_ref(pixel, angle, m_xmt_aperture.data->m_radius,
																 xmt_apo.data->m_ref);
								}
								else {
									xmt_apo.data->m_ref[0] = recv_offset;
									xmt_apo.data->m_ref[1] = 0.0;
									xmt_apo.data->m_ref[2] = 0.0;
								}
								/* Added because of introduction of m_apo_dir for
									 transmit apodization */
								xmt_apo.data->m_apo_dir[0] = -sin(angle);
								xmt_apo.data->m_apo_dir[1] = 0;
								xmt_apo.data->m_apo_dir[2] = cos(angle);

                float_type val_rcv = float_type(1.0);
                float_type val_xmt = float_type(1.0);           

                if (xmt_apo.fixed()) {
									val_xmt = xmt_apo.fixApo(trans_no, pixel, xmt_ap_size);
									// IEEE
									//val_xmt = rcv_apo.fixApo(5+trans_no, pixel, xmt_ap_size);
                }
                if (xmt_apo.dynamic()) {
                  xmt_apo.setPosition(pixel);
                  xmt_apo.setFnumber();
                  val_xmt = val_xmt * xmt_apo.dynApo(trans_no,pixel);
                }
                result = result*val_xmt;                
          
                if (rcv_apo.fixed()) {
                  val_rcv = rcv_apo.fixApo(i_rcv, pixel, rcv_ap_size);
                }
                if (xmt_apo.dynamic()) {
                  rcv_apo.setPosition(pixel);
                  rcv_apo.setFnumber();
                  val_rcv = val_rcv * rcv_apo.dynApo(i_rcv, pixel);
                }
                result = result*val_rcv;
                image[i_x*nz + i_z]  += result; // Add a lot of zeros
              } // end for i_angle
            } // end if findex < (n_rf_samples-1)
          } // end for i_rcv
        } // end for i_x
      } // end for i_y
    } // end for i_z
    for (i=0; i < n_rcv_elements ; i++) {
      delete interpl[i];
    }
    delete[] interpl;
  } // end for i_xmt
#if defined(HAVE_PTHREAD_H)
  pthread_exit(NULL);
#else
  return 0;
#endif
}

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
