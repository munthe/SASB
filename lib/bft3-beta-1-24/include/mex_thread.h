/***************************************************************
 *                                                             *
 * Project              SIMD BEAMFORM MEX                      *
 * Module               THREAD UTILITY FUNCTIONS               *
 *                                                             *
 * $Id: mex_thread.h,v 1.21 2011-07-27 22:56:02 jmh Exp $ 
 *                                                             *
 * $Author: jmh $                                              *
 *                                                             *
 * $Date: 2011-07-27 22:56:02 $                                *
 *                                                             *
 * $State: Exp $                                               *
 *                                                             *
 *-------------------------------------------------------------*
 * $Log: mex_thread.h,v $
 * Revision 1.21  2011-07-27 22:56:02  jmh
 * *** empty log message ***
 *
 * Revision 1.20  2011/04/17 11:10:23  jmh
 * *** empty log message ***
 *
 * Revision 1.19  2011/04/03 23:23:06  jmh
 * *** empty log message ***
 *
 * Revision 1.18  2011/02/27 22:24:36  jmh
 * *** empty log message ***
 *
 * Revision 1.17  2010/05/02 13:38:39  jmh
 * *** empty log message ***
 *
 * Revision 1.16  2010/04/04 13:22:24  jmh
 * Experiment with signal handler
 *
 * Revision 1.15  2010/04/04 11:46:52  jmh
 * *** empty log message ***
 *
 * Revision 1.14  2010/04/01 01:32:50  jmh
 * *** empty log message ***
 *
 * Revision 1.13  2010/03/31 21:31:33  jmh
 * *** empty log message ***
 *
 * Revision 1.12  2010/03/30 19:34:15  jmh
 * VS XMT and RCV - slow but it works
 *
 * Revision 1.11  2010/03/27 19:54:12  jmh
 * Main thread initializes interpolators
 *
 * Revision 1.10  2010/03/26 15:39:53  jmh
 * Experiments with pthreads on blows
 *
 * Revision 1.9  2010/03/23 14:50:34  jmh
 * *** empty log message ***
 *
 * Revision 1.8  2010/03/19 17:52:26  jmh
 * AV when beamforming t1.m
 *
 * Revision 1.7  2010/03/18 01:04:44  jmh
 * doxygen
 *
 * Revision 1.6  2010/03/18 00:49:51  jmh
 * mex_thread.h ready for use, TODO: Signal handling
 *
 * Revision 1.5  2010/03/18 00:26:50  jmh
 * Solved threading issues of members
 *
 * Revision 1.4  2010/03/14 11:53:34  jmh
 * Pretty print
 *
 * Revision 1.3  2010/03/12 23:33:49  jmh
 * Added 32-bit linux Makefile
 *
 * Revision 1.2  2010/03/12 15:48:53  jmh
 * *** empty log message ***
 *
 * Revision 1.1  2010/02/09 21:25:51  jmh
 * *** empty log message ***
 *
 ***************************************************************/

/***************************************************************
 * TODO:
 *  - Integrate event handler
 *  - Unit test before integration with BFT
 *  -
 ***************************************************************/

#ifndef _MEX_THREAD_H
#define _MEX_THREAD_H

// Assume pthreads available
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef _WIN32
# ifndef NOMINMAX
#  define NOMINMAX
# endif
# include <windows.h>
# include <process.h>
# include <io.h>
# if defined(HAVE_PTHREAD_H)
#  include <pthread.h>
# endif
#else
# if defined(HAVE_SIGNAL_H)
#  include <csignal>
# endif
# if defined(HAVE_PTHREAD_H)
#  include <pthread.h>
# endif
# include <unistd.h>
# include <sys/syscall.h>
#endif

#if defined(HAVE_PTHREAD_H)
template<class T, void*(T::*thread_func)(void*)>
#elif defined(_WIN32)
template<class T, unsigned int(__stdcall T::*thread_func)(void*)>
#endif

class pthread_launcher {
 public:
  pthread_launcher(T* obj=NULL, void* arg=NULL) : _obj(obj), _arg(arg) {}
#if defined(HAVE_PTHREAD_H)
  void *launch() { return (_obj->*thread_func)(_arg);}
#elif defined(_WIN32)
  unsigned int launch() { return (_obj->*thread_func)(_arg);}
#endif
 private:
  /// Object pointer
  T* _obj;
  /// Command argument for member function
  void *_arg;
};

// Launch thread function
template<class T>
#if defined(HAVE_PTHREAD_H)
void *launch_member_function(void *obj)
#elif defined(_WIN32)
unsigned int __stdcall launch_member_function(void *obj)
#endif
{
  T* launcher = reinterpret_cast<T*>(obj);
  return launcher->launch();
}

#endif // _MEX_THREAD_H
