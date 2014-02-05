/******************************************************************************
 *                                                                            *
 * Project              SIMD BEAMFORM MEX                                     *
 * Module               Image Class                                           *
 *                                                                            *
 * $Id: image.cpp,v 1.178 2011-08-02 13:44:07 jmh Exp $
 *                                                                            *
 * $Author: jmh $                                                             *
 *                                                                            *
 * $Date: 2011-08-02 13:44:07 $                                               *
 *                                                                            *
 * $State: Exp $                                                              *
 *----------------------------------------------------------------------------*
 *                                                                            *
 *****************************************************************************/

/****************************************************************************
 * TODO:
 *  - Introduce messages
 *  - Make orientation for element on emit_aperture and apodize in transmit
      according to distance to emission and distance to this line
 *  - Remove ppwave for receive aperture
 ****************************************************************************/

/*
  A call to the signal() function defines a signal handler function to
  be given control in the event a particular signal is raised. The C
  standard constrains the behavior of this handler function in a number
  of respects -- it can only terminate in certain ways, it usually can't
  call other standard library functions, it can't refer to static
  objects unless they are of type volatile sig_atomic_t, and so forth
  (see Section 7.7.1.1 of the ISO C standard for full details). Any
  violation of these constraints results in undefined behavior.

  This left open the question of what constraints a C++ signal handling
  function has. The standards committee resolved this during the past
  year by adding language to the standard defining the concept of a
  "plain old function" (POF), analogous to the existing concept of plain
  old data (POD). A POF is a function that only uses the common subset
  of the C and C++ languages. A C++ signal handler only has defined
  behavior if it is a POF and if it would have defined behavior under
  the C standard; in particular, any handler which is not a POF -- i.e.,
  which uses any C++ features -- will have undefined behavior.

*/

// Debugging
#include <cstdio>
#include <cstdlib>

// Error handling
#include <cerrno>
#include <cstring> // strerror

#include "aperture.h"
#include "apodization.h" // screws up header hierachy
#include "image.h"
#include "sample_interpolate.h"

#ifdef HAVE_MQUEUE_H
#undef HAVE_MQUEUE_H
#endif

#ifndef MSGMAX
# define MSGMAX 8192
#endif

# ifdef __cplusplus 
extern "C" bool utIsInterruptPending();
extern "C" void utSetInterruptPending(bool);
# else
extern bool utIsInterruptPending();
extern void utSetInterruptPending(bool);
# endif

#ifdef __linux__
static pid_t gettid( void ) {
  pid_t pid;
  CallErr(pid = syscall, (__NR_gettid));
  return pid;
}
#endif

// Template instantiation
#ifdef HAVE_MQUEUE_H
template class pthread_launcher<Image,&Image::thread_func>;
#else
template class pthread_launcher<Image,&Image::beamform_threaded>;
#endif

// Static variables
size_t Image::nthreads = 1;
size_t Image::DataRep::ref_cnt = 0; // TODO: Use for stopping threads
																		// when last object is deleted
Image::thread_arg_t Image::threadarg[N_MAX_THREADS];
#ifdef HAVE_PTHREAD_H
pthread_t Image::threads[N_MAX_THREADS];
pthread_attr_t Image::attr = {{0}};
#endif

#ifdef HAVE_MQUEUE_H
// Static variables
bool  Image::threads_initialized = false;
mqd_t Image::mqd_master = 0;
mqd_t Image::mqd_client = 0;

int Image::mq_clear(const char* qname) {
  mqd_t mqd;
  struct mq_attr qattr, old_qattr;
  unsigned int prio = 0;
  char buf[MSGMAX];

  // TODO: Change O_RDWR | O_CREAT to O_RDONLY
  CallErrExit(mqd = mq_open,
              (qname, O_RDWR | O_CREAT,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
               NULL), EXIT_SUCCESS);
  
  // Remove messages on the queue
  mq_getattr(mqd, &qattr);

  if (qattr.mq_curmsgs !=0) {
    qattr.mq_flags = O_NONBLOCK;
    mq_setattr (mqd, &qattr, &old_qattr);    
    while (mq_receive(mqd, &buf[0], qattr.mq_msgsize, &prio) != -1) {
    }
    if (errno != EAGAIN) { 
      perror ("mq_receive()");
      return EXIT_FAILURE;
    }
    // Restore attributes
    mq_setattr(mqd, &old_qattr, 0);            
  }
  CallErrExit(mq_close,(mqd),EXIT_FAILURE);
  return EXIT_SUCCESS;
}
#endif

Image::Image() : data(NULL) {
}

Image::Image(const Image &other) {
  data = other.data;
  if (data)
    data->m_cnt++;
}

Image& Image::operator=(const Image &other) {
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

Image::~Image() {
  if (data)
    if (--data->m_cnt == 0) {
      delete data;
		}
}

Image::Image(const Aperture& Th_t, const Aperture& Th_r,
             const vector<Apodization>& Ah_t,
             const vector<Apodization>& Ah_r,
             const vector<Line>& lines) : data(NULL) {
	
	data = new Image::DataRep(Th_t, Th_r, Ah_t, Ah_r, lines);
  data->m_interp_type = Image::linear;
}

#if (defined(HAVE_PTHREAD_H) && defined(HAVE_MQUEUE_H))
void* Image::thread_func(void *ptarg) {

  char buf[MSGMAX];
  thread_arg_t* threadarg = NULL;
  cpu_set_t set;
  struct mq_attr qattr;
  unsigned int prio;

  threadarg = (thread_arg_t*) ptarg;

  int thread_id = threadarg->thread_id;
	
	// Priority set to thread ID
  prio = thread_id;

#ifdef __linux__
  CPU_ZERO( &set );
  CPU_SET( threadarg->cpu_id, &set );

	/* Do not explicitly set processor id
  CallErrExit(sched_setaffinity,
              (gettid(), sizeof(cpu_set_t), &set),NULL);
	*/
#elif defined(_WIN32)
	HANDLE hThread = GetCurrentThread(void);
	SetThreadIdealProcessor(hThread,threadarg->cpu_id);
#endif

  mqd_t lmqd_master;
  mqd_t lmqd_client;
  CallErrExit(lmqd_master = mq_open,
              ("/jmh-master",O_WRONLY,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
               NULL), NULL);

  // Post Message
  if (mq_send(lmqd_master,"READY",5,prio)==-1)
    perror ("mq_send(): READY");

  int state = 0; // Ready
  void* thread_retval = NULL;
  

  CallErrExit(lmqd_client = mq_open,
              ("/jmh-client",O_RDONLY,
               S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
               NULL), NULL);

  mq_getattr (lmqd_client, &qattr);

  qattr.mq_flags &= ~O_NONBLOCK;
  mq_setattr(lmqd_client,&qattr,0);
  int rett = 0;

  while(mq_receive (lmqd_client, &buf[0], qattr.mq_msgsize, &prio) != -1) {
    if (errno == EAGAIN) { 
      printf("Apparently non-blocking\n");
      pthread_exit(NULL);
    }

    if (!strcmp(buf,"RESET"))
      state = 0;
    else if (!strncmp(buf,"RUN",3)) {
      state = 1;
    }
    else if (!strncmp(buf,"EXIT",4)) {
      state = 3;
      break;
    }
    else {
      printf("Unknown message: %s", buf);
    }

    switch (state) {
    case 0:
    case 1:
      thread_retval = beamform_threaded((void*)threadarg);
      state = 2;
    case 2:
      rett = mq_send(lmqd_master,"DONE",4,prio);
      if (rett==-1) {
        printf("Error sending DONE\n");
        perror ("mq_send()");
      }
      state = 0;
      break; // break-out of switch only
    }
  };

	// TODO: Make CallErrExit call pthread_exit(NULL)
  CallErrExit(mq_close,(lmqd_client),NULL);
  CallErrExit(mq_close,(lmqd_master),NULL);
  pthread_exit(NULL);
}
#endif

#if defined(HAVE_PTHREAD_H)
void* Image::beamform_threaded(void *ptarg)
#else
  unsigned __stdcall Image::beamform_threaded(void *ptarg)
#endif
{

  ///////////////////////////////////////////////////
  // Local variables
  ///////////////////////////////////////////////////
  const thread_arg_t *arg = (thread_arg_t*) ptarg;
  const size_t n_rf_samples   = arg->n_rf_samples;
  const size_t start_line     = arg->start_line;
  const size_t end_line       = arg->end_line;
  const float_type* rf_data   = arg->rf_data;
  const size_t n_rcv_elements = data->m_rcv_aperture.n_elements();

  // Receive delays (constant)
  const float_type* rcv_delays = data->m_rcv_aperture.data->m_delays;

  // Virtual sources
  const float_type* xmt_focus  = data->m_xmt_aperture.data->m_focus;
  const float_type* rcv_focus  = data->m_rcv_aperture.data->m_focus;

  // System variables (global)
  const float_type c          = *data->m_xmt_aperture.data->m_c;
  const float_type fs         = *data->m_xmt_aperture.data->m_fs;

  // Output image
  float_type* image = arg->image;

  ///////////////////////////////////////////////////
  // Temporary variables
  ///////////////////////////////////////////////////
  size_t trans_no;
  float_type pixel_x0, pixel_y0, pixel_z0;

  // Pixel location
  ALIGN16_BEGIN float_type pixel[3] ALIGN16_END;

  // Rcv element position
  ALIGN16_BEGIN float_type xyz_r[3] ALIGN16_END;

  // Virtual sources
  ALIGN16_BEGIN float_type xyz_rf[3] ALIGN16_END;
  ALIGN16_BEGIN float_type xyz_tf[3] ALIGN16_END;

  // Normal vectors used for pp-wave support
  ALIGN16_BEGIN float_type xyz_tn[3] ALIGN16_END;
  ALIGN16_BEGIN float_type xyz_rn[3] ALIGN16_END;

  ALIGN16_BEGIN float_type euler_t[3] ALIGN16_END;

  float_type dx, dy, dz;
  float_type dist_pixel_to_xmt;
  float_type dist_pixel_to_rcv;

  // Distance from xmt focus (VS) to xmt reference
  float_type dist_focus_to_xmt;

  // Distance from rcv focus (VS) to rcv element
  float_type dist_focus_to_rcv;

  // Floating point sample index
  float_type findex;

  memset(xyz_tn,0,3*sizeof(float_type));
  memset(xyz_rn,0,3*sizeof(float_type));

#ifdef __linux__
# ifndef HAVE_MQUEUE_H
  cpu_set_t set;
  CPU_ZERO(&set);
  CPU_SET(arg->cpu_id, &set);
  CallErrExit(sched_setaffinity,
              (gettid(), sizeof(cpu_set_t), &set ),
              NULL);
# endif
#elif defined(_WIN32)
	HANDLE hThread = GetCurrentThread(void);
	SetThreadIdealProcessor(hThread,threadarg->cpu_id);
#endif

  ALIGN16_BEGIN float_type xyz_t[3] ALIGN16_END;
  SampleInterpolate<size_t,float_type>** interpl = NULL;

  // TODO: Use inter-thread derivatives for spline interpolation
  for (size_t i_emission=0 ; 
			 i_emission<data->m_xmt_aperture.data->m_nemissions ; i_emission++) {

    interpl =
      new SampleInterpolate<size_t,float_type>*[n_rcv_elements];
    for (size_t i_rcv=0; i_rcv < n_rcv_elements ; i_rcv++) {

      interpl[i_rcv] =
        new SampleInterpolate<size_t, float_type>
        (&rf_data[i_rcv*n_rf_samples+i_emission*n_rf_samples*n_rcv_elements],
         n_rf_samples);
      interpl[i_rcv]->setMethod(data->m_interp_type);
    }
    xyz_t[0] = data->m_xmt_aperture.data->m_center_focus[0+i_emission*3];
    xyz_t[1] = data->m_xmt_aperture.data->m_center_focus[1+i_emission*3];
    xyz_t[2] = data->m_xmt_aperture.data->m_center_focus[2+i_emission*3];

    euler_t[0] = data->m_xmt_aperture.data->m_euler[0+i_emission*3];
    euler_t[1] = data->m_xmt_aperture.data->m_euler[1+i_emission*3];
    euler_t[2] = data->m_xmt_aperture.data->m_euler[2+i_emission*3];

    trans_no = (arg->trans_no[i_emission]-1);

    if (xmt_focus) {
      // Virtual source
      xyz_tf[0] = xmt_focus[0];
      xyz_tf[1] = xmt_focus[1];
      xyz_tf[2] = xmt_focus[2];
    }
    else {
      // Virtual source equals center_focus
      xyz_tf[0] = xyz_t[0];
      xyz_tf[1] = xyz_t[1];
      xyz_tf[2] = xyz_t[2];
    }
    
    if (rcv_focus) {
      // Virtual source
      xyz_rf[0] = rcv_focus[0];
      xyz_rf[1] = rcv_focus[1];
      xyz_rf[2] = rcv_focus[2];
    }
    else {
			// Virtual source (not used)
      xyz_rf[0] = float_type(0);
      xyz_rf[1] = float_type(0);
      xyz_rf[2] = float_type(0);
    }
    
    for (size_t i_line=start_line ; i_line < end_line ; i_line++) {
      
      Line line = data->m_lines[i_line];
      Apodization rcv_apodization = data->m_rcv_apo[i_line];
      Apodization xmt_apodization = data->m_xmt_apo[i_line];

      // TODO: Specialize for pp-wave
      if (data->m_xmt_aperture.ppwave()) {
        basis_vector(xyz_tn,euler_t,2);
      }
      if (data->m_rcv_aperture.ppwave()) {
        basis_vector(xyz_rn,rcv_apodization.data->m_euler,2);
      }
      
      // Using delay to skip beginning of lines
      line.findMinSample(arg->delay[i_emission],data->m_xmt_aperture);
      
      line.setApodization2(xmt_apodization,
                           rcv_apodization, trans_no);
      
      dist_focus_to_xmt =
        sqrt(SQUARE(xyz_t[0]-xyz_tf[0]) +
             SQUARE(xyz_t[1]-xyz_tf[1]) +
             SQUARE(xyz_t[2]-xyz_tf[2]));
      dist_focus_to_xmt = fs*dist_focus_to_xmt/c;

      pixel_x0 = line.origin(0);
      pixel_y0 = line.origin(1);
      pixel_z0 = line.origin(2);
      dx = line.dx();
      dy = line.dy();
      dz = line.dz();

      size_t line_length = (size_t) floor(line.length() / line.dr());

      // Multiple focal points
      for (size_t i_limit=0; (i_limit+1) < line.data->m_limits.size() ;
					 i_limit++) {

        size_t line_start = line.data->m_limits[i_limit];
        size_t line_end = line.data->m_limits[i_limit+1];

        for (size_t i_sample=line_start ; i_sample < line_end ; i_sample++) {

          pixel[0] = pixel_x0 + i_sample*dx;
          pixel[1] = pixel_y0 + i_sample*dy;
          pixel[2] = pixel_z0 + i_sample*dz;

          // Critical (TODO: Template specialize function)
          if (data->m_xmt_aperture.ppwave()) {
            dist_pixel_to_xmt =
              fabs((pixel[0] - xyz_tf[0])*xyz_tn[0] +
                   (pixel[1] - xyz_tf[1])*xyz_tn[1] +
                   (pixel[2] - xyz_tf[2])*xyz_tn[2]);
          }
          else {
            dist_pixel_to_xmt =
              sqrt(SQUARE(xyz_tf[0]-pixel[0]) +
                   SQUARE(xyz_tf[1]-pixel[1]) +
                   SQUARE(xyz_tf[2]-pixel[2]));
          }

          // Critical (TODO: Template specialize function)
          if (rcv_apodization.dynamic()) {
            rcv_apodization.setPosition(pixel);
            rcv_apodization.setFnumber();
          }

          // Critical (TODO: Template specialize function)
          if (xmt_apodization.dynamic()) {
            xmt_apodization.setPosition(pixel);
            xmt_apodization.setFnumber();
          }

          // Beamformation
          dist_pixel_to_xmt = fs*dist_pixel_to_xmt/c;

          if (xmt_focus) {
            if (line.data->m_xmt_vs_order[i_limit] == 0)
              dist_pixel_to_xmt = dist_focus_to_xmt - dist_pixel_to_xmt;
            else
              dist_pixel_to_xmt = dist_focus_to_xmt + dist_pixel_to_xmt;
          }
          dist_pixel_to_xmt -= arg->delay[i_emission]*fs;

          // TODO: Test with this as the outermost loop
          for (size_t i_rcv=0 ; i_rcv < n_rcv_elements ; i_rcv++) {

            // Distance to receive
            xyz_r[0] = data->m_rcv_aperture(0,i_rcv);
            xyz_r[1] = data->m_rcv_aperture(1,i_rcv);
            xyz_r[2] = data->m_rcv_aperture(2,i_rcv);

            // Critical (TODO: Specialize function)
            if (rcv_focus) {
              xyz_rf[0] = rcv_focus[0];
              xyz_rf[1] = rcv_focus[1];
              xyz_rf[2] = rcv_focus[2];          
            }
            else {
              xyz_rf[0] = xyz_r[0];
              xyz_rf[1] = xyz_r[1];
              xyz_rf[2] = xyz_r[2];
            }

            // focus_rcv -> ref_rcv (when no rcv_vs, this is zero)
            dist_focus_to_rcv =
              sqrt(SQUARE(xyz_r[0]-xyz_rf[0]) +
                   SQUARE(xyz_r[1]-xyz_rf[1]) +
                   SQUARE(xyz_r[2]-xyz_rf[2]));

            if (data->m_rcv_aperture.ppwave()) {
              dist_pixel_to_rcv =
                fabs((pixel[0] - xyz_rf[0])*xyz_rn[0] +
                     (pixel[1] - xyz_rf[1])*xyz_rn[1] +
                     (pixel[2] - xyz_rf[2])*xyz_rn[2]);
            }
            else {
              dist_pixel_to_rcv =
                sqrt(SQUARE(xyz_rf[0]-pixel[0]) +
                     SQUARE(xyz_rf[1]-pixel[1]) +
                     SQUARE(xyz_rf[2]-pixel[2]));
            }

            if (rcv_focus) {
              if (line.data->m_rcv_vs_order[i_limit] == 0)
                dist_pixel_to_rcv =
                  dist_focus_to_rcv - dist_pixel_to_rcv;
              else
                dist_pixel_to_rcv =
                  dist_focus_to_rcv + dist_pixel_to_rcv;
            }
            dist_pixel_to_rcv = fs*dist_pixel_to_rcv/c;

            findex = dist_pixel_to_rcv + dist_pixel_to_xmt;

            // Constant delays
            if (rcv_delays) {
              findex = findex + rcv_delays[i_rcv]*fs;
            }

            if ((findex > 0.0) && (findex < (n_rf_samples-1))) {

              float_type val = float_type(0.0);     
              float_type val2 = float_type(0.0);    

              float_type result = float_type(0.0);

              result = (*interpl[i_rcv])(findex);

              // Critical (TODO: Specialize function)
              if (xmt_apodization.dynamic()) {
                val2 = xmt_apodization.dynApo(trans_no,pixel);      
                result = result * val2;
              }
              if (rcv_apodization.dynamic()) {
                val = rcv_apodization.dynApo(i_rcv,pixel);
                result = result * val;
              }
              if (xmt_apodization.manual()) {
                result = result *
                  xmt_apodization(trans_no,line.data->m_xmt_apo_order[i_limit]);
              }
              if (rcv_apodization.manual()) {
                result = result *
                  rcv_apodization(i_rcv,line.data->m_rcv_apo_order[i_limit]);
              }
              image[i_line*line_length + i_sample]  += result;
            }
          } // Receive elements
        } // Samples along line
      } // Multiple focus points
      if (utIsInterruptPending()) { /* check for a Ctrl-C event */
				utSetInterruptPending(false);
        for (size_t i_rcv=0; i_rcv < n_rcv_elements ; i_rcv++) {
          delete interpl[i_rcv];
        }
        delete[] interpl;
# if defined(HAVE_PTHREAD_H)
#  if (HAVE_MQUEUE_H)
				return NULL;
#  else
        pthread_exit(NULL);
#  endif
# elif defined(_WIN32)
        return 0;
# endif
      }
    } // Lines
    // Clean-up interpolation objects
    for (size_t i_rcv=0; i_rcv < n_rcv_elements ; i_rcv++) {
      delete interpl[i_rcv];
    }
    delete[] interpl;
  } // Emissions

#if defined(HAVE_PTHREAD_H)
#  if (HAVE_MQUEUE_H)
				return NULL;
#  else
        pthread_exit(NULL);
#  endif
#elif defined(_WIN32)
  return 0;
#endif
}

// RF dimensions (rf_samples, rcv_channels, xmt_channels)
bool Image::Beamform(float_type* image, const float_type* rf_data,
                     const float_type* delay, const size_t n_rf_samples,
                     const uint32_t* trans_no) {

  size_t n_lines = data->m_lines.size();
  bool retval = true;

#ifndef HAVE_THREAD
  // Single threaded (TODO: Verify)
  thread_arg_t threadarg;
  threadarg.start_line = 0;
  threadarg.end_line = n_lines;
  threadarg.image = image;
  threadarg.rf_data = rf_data;
  threadarg.delay = delay;
  threadarg.n_rf_samples = n_rf_samples;
  threadarg.trans_no = trans_no;
  threadarg.thread_id = 0;
  threadarg.emission_no = 0;
  void* thread_retval = beamform_threaded((void*)&threadarg);
#else

  size_t i;
  int nproc;

# if defined(HAVE_PTHREAD_H)
# elif defined(_WIN32)
  unsigned int threadID;
  HANDLE threads[N_MAX_THREADS];
# endif

# if defined(__linux__)
  cpu_set_t cpuset;
  CPU_ZERO( &cpuset );
  CallErr(sched_getaffinity,
          (gettid(), sizeof( cpu_set_t ), &cpuset));

  nproc = CPU_COUNT(&cpuset);
# elif defined(_WIN32)
  SYSTEM_INFO info;
  GetSystemInfo(&info);
  nproc = info.dwNumberOfProcessors;
	DWORD sysMask;
	DWORD procMask;
	HANDLE hProc = GetCurrentProcess(void);
  GetProcessAffinityMask(hProc,&procMask,&sysMask);
# endif

# ifdef HAVE_MQUEUE_H
  pthread_launcher<Image,
    &Image::thread_func> launcher[N_MAX_THREADS];
# else
  pthread_launcher<Image,
    &Image::beamform_threaded> launcher[N_MAX_THREADS];
# endif

  setbuf(stdout, NULL);
  setbuf(stderr, NULL);

# ifdef HAVE_PTHREAD_H
  CallErr(pthread_attr_init,(&Image::attr));
#  ifdef HAVE_MQUEUE_H
  CallErr(pthread_attr_setdetachstate, (&Image::attr, PTHREAD_CREATE_DETACHED));
#  else
  CallErr(pthread_attr_setdetachstate, (&Image::attr, PTHREAD_CREATE_JOINABLE));
#  endif
# endif

	// Populate structs for threads
  for (i=0;i<nthreads;i++) {
    threadarg[i].start_line   = 0+i*(n_lines/nthreads);
    threadarg[i].end_line     = (n_lines/nthreads)+i*(n_lines/nthreads);
    threadarg[i].image        = image;
    threadarg[i].rf_data      = rf_data;
    threadarg[i].delay        = delay;
    threadarg[i].n_rf_samples = n_rf_samples;
    threadarg[i].trans_no     = trans_no;
    threadarg[i].thread_id    = i;
    threadarg[i].cpu_id       = ((int) i) % nproc;
		if (i==(nthreads-1))
			threadarg[i].end_line   = n_lines;      
  }

# ifdef HAVE_MQUEUE_H
  // Mesage queues
  struct mq_attr qattr;
  
  char buf[MSGMAX];
  unsigned int prio;

	// Create two message queues
	CallErrExit(Image::mqd_master = mq_open,
							("/jmh-master", O_RDWR | O_CREAT,
							 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
							 NULL), EXIT_FAILURE);
	CallErrExit(mq_close,(Image::mqd_master),EXIT_FAILURE);

	CallErrExit(Image::mqd_client = mq_open,
							("/jmh-master", O_RDWR | O_CREAT,
							 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
							 NULL), EXIT_FAILURE);
	CallErrExit(mq_close,(Image::mqd_client),EXIT_FAILURE);

	/* Initiate launcher, note that the threadarg are global and updated
		 in between calls to beamform_threaded */
  for (i=0;i<nthreads;i++) {
    launcher[i]               =
      pthread_launcher<Image,&Image::thread_func>(this,
                                                  &threadarg[i]);
  }

  if ((!Image::threads_initialized) && !(nthreads>N_MAX_THREADS)) {

    /* Clear message queues (and create if they don't exist), TODO: Do
			 this in two stages */
    CallErrExit(mq_clear,("/jmh-master"),false);
    CallErrExit(mq_clear,("/jmh-client"),false);

		/* Initiate threads */
    for (i=0;i<nthreads;i++) {
      CallErr(pthread_create,
              (&threads[i], &Image::attr,
               launch_member_function<pthread_launcher<Image,
               &Image::thread_func> >, 
               &launcher[i]));
    }
    Image::threads_initialized = true;
	}

	/* Signal idling threads to run */
	CallErrExit(Image::mqd_client = mq_open,
							("/jmh-client", O_WRONLY,
							 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
							 NULL), false);
	for (i=0;i<nthreads;i++) {
		prio = i;
		if (mq_send(mqd_client,"RUN",3,prio)==-1)
			perror ("mq_send()");
	}
	mq_close(Image::mqd_client);
	
	/* Wait for threads to finish */
	CallErrExit(Image::mqd_master = mq_open,
							("/jmh-master", O_RDONLY,
							 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
							 NULL), false);
	
	mq_getattr(Image::mqd_master, &qattr);
	
	size_t n_to_go = nthreads;
	size_t n_idle  = 0;
	/* Wait for threads to finish */
	while (true) {
		mq_receive(Image::mqd_master, &buf[0], qattr.mq_msgsize, &prio);
		if (!strncmp(buf,"DONE",4)) {
			n_to_go--;
		}
		else if (!strncmp(buf,"READY",5)) {
			n_idle++;
		}
		else {
			printf("Unexpected message: %s\n",buf);
			retval = false;
			break;
		}
		if (n_to_go==0)
			break;
	}
	mq_close(Image::mqd_master);
# else
  /* Without message queues (slower) */
  for (i=0;i<nthreads;i++) {
    launcher[i]               =
      pthread_launcher<Image,&Image::beamform_threaded>(this,
                                                        &threadarg[i]);
  }
  for (i=0;i<nthreads;i++) {
#  if defined(HAVE_PTHREAD_H)
    CallErr(pthread_create,
            (&threads[i], &Image::attr,
             launch_member_function<pthread_launcher<Image,
             &Image::beamform_threaded> >, 
             &launcher[i]));
#  elif defined(_WIN32)
    threads[i] =
      (HANDLE)_beginthreadex(NULL, 0,
                             launch_member_function<pthread_launcher<Image,
                             &Image::beamform_threaded> >,
                             &launcher[i], 0, &threadID );
#  endif
  }
  
  for (i = 0; i < nthreads; i++) {
#  if defined(HAVE_PTHREAD_H)
    CallErr(pthread_join,(threads[i],NULL));
#  elif defined(_WIN32)
    WaitForSingleObject(threads[i], INFINITE );
#  endif
  }

#  if defined(HAVE_PTHREAD_H)
  CallErr(pthread_attr_destroy,(&Image::attr));
#  endif
# endif
#endif
  return retval;
}

Image::DataRep::DataRep(const Aperture& xmt_aperture, 
												const Aperture& rcv_aperture,
												const vector<Apodization>& xmt_apo,
												const vector<Apodization>& rcv_apo,
												const vector<Line>& lines) : m_cnt(1) {
	m_xmt_aperture = xmt_aperture;
	m_rcv_aperture = rcv_aperture;
	m_xmt_apo = xmt_apo;
	m_rcv_apo = rcv_apo;
	m_lines = lines;
}

// Destructor
Image::DataRep::~DataRep() {
#ifdef HAVE_MQUEUE_H
	size_t i;
  unsigned int prio = 0;

	if (Image::threads_initialized) {
		// Signal threads to exit
		CallErr(Image::mqd_client = mq_open,
						("/jmh-client", O_WRONLY | O_CREAT,
						 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH,
						 NULL));
	
	for (i=0;i<nthreads;i++) {
		prio = i;
		if (mq_send(mqd_client,"EXIT",4,prio)==-1)
			perror ("mq_send()");
	}
	
	// Wait for threads to exit
	for (i = 0; i < nthreads; i++) {
		CallErr(pthread_join,(threads[i],NULL));
	}
	mq_close(Image::mqd_client);
	
	Image::threads_initialized = false;
# if defined(HAVE_PTHREAD_H)
	CallErr(pthread_attr_destroy,(&attr));
# endif
	}
#endif
}

/* Local variables: */
/* default-tab-width: 2 */
/* End: */
