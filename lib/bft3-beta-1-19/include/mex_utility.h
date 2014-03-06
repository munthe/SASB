/*************************************************************** 
 *                                                             * 
 * Project              SIMD BEAMFORM MEX                      * 
 * Module               AUTOPTR AND GARBAGE COLLECTOR          * 
 *                                                             * 
 * $Id: mex_utility.h,v 1.33 2011-04-17 12:49:20 jmh Exp $   
 *                                                             * 
 * $Author: jmh $                                              * 
 *                                                             * 
 * $Date: 2011-04-17 12:49:20 $                                * 
 *                                                             * 
 * $State: Exp $                                               * 
 *                                                             * 
 * Adapted from a design by Mike Stevens.                      * 
 *                                                             * 
 *-------------------------------------------------------------*
 *                                                             *
 *
 ***************************************************************/

/***************************************************************
 * TODO:
 *  - Branch on sprintf(size_t)
 *  -
 ***************************************************************/

#ifndef MEX_AUTO_PTR_H
#define MEX_AUTO_PTR_H

#include <cstdio>
#include <cstdlib>

#include <typeinfo>
#include <list>
#include <cstring>

#ifdef _MSC_VER
 #pragma warning( disable : 4996 ) // suppress warning: 'sprintf' was declared deprecated
 #define _CRT_SECURE_NO_WARNINGS
#endif

#if defined (_WIN32)
# ifndef NOMINMAX
#  define NOMINMAX
# endif
# include <windows.h>
#elif defined (__linux__)
# define __STDC_FORMAT_MACROS 1
# include <inttypes.h>
# include <unistd.h>
extern const char *__progname;
#endif

#include <mex.h>

#ifndef mxPOINTER_CLASS
# define mxPOINTER_CLASS mxUINT64_CLASS
#endif

#ifndef mxIsPointer
# define mxIsPointer mxIsScalarUInt64
#endif

#ifndef mxPOINTER_TYPE
# define mxPOINTER_TYPE UINT64_T
#endif

#ifndef mxIsPointerArray
# define mxIsPointerArray mxIsUint64
#endif

/*///////////////////////////////////////////////////////////////////
//
// Name:
//   MexAutoPtr
//
// Decription:
//
//   Object wrapper class similar to STL auto_ptr
//
///////////////////////////////////////////////////////////////////*/

#if (defined(__GNUC__) && __GNUC__ >= 3) \
 || (defined(_WIN32)) \
 || (defined(linux) && defined(__INTEL_COMPILER) && defined(__ICC))
# define MEX_TYPE_ID_NAME
#endif 

// For this compiler at least, cross-shared-library type_info
// comparisons don't work, so use typeid(x).name() instead. It's not
// yet clear what the best default strategy is.
#ifdef MEX_TYPE_ID_NAME
  typedef char const* base_id_t;
#else
  typedef std::type_info const* base_id_t;
#endif
 
// Forward declaration of Garbage Collector
template<typename T> class GarbageCollector;

template <typename T>
class MexAutoPtr {
 public:
  // Auto pointer ctor for free-store allocated objects. The auto
  // pointer takes ownership, and will delete object when it is
  // destroyed. Be sure not to delete any aliases to the object
  // argument
#ifdef MEX_TYPE_ID_NAME
  MexAutoPtr(T*& ptr) : type(typeid(T).name()), t(ptr) { 
#else
  MexAutoPtr(T*& ptr) : type(&typeid(T)), t(ptr)  { 
#endif
    signature = this; 
    GarbageCollector<T>::register_handle(this);
    ptr = 0;
  } 

  ~MexAutoPtr() { 
    delete t;        // destroy object
    signature= NULL; // destroy signature
  } 

  // Convert mxArray (passed to mex-function) to an MexAutoPtr<T>.
  static MexAutoPtr* from_mex_handle( const mxArray* ma );

  // Convert MexAutoPtr<T> to a mxArray handle (to pass back from mex-function).
  mxArray* to_mex_handle(); 

  // Get the actual object contained by handle
  T& get_object() const { return *t; }

 private:

  // Temporarily public
  public:
  MexAutoPtr* signature;      // used as a unique object signature 

  base_id_t type; // type checkig information
  T *t; // object pointer

  friend class GarbageCollector<T>;
};

// Helper functions for simple usage

// Create mex handle to object t (where t is heap allocated). 
// Client no longer owns t, and so must not delete it.
template <typename T> mxArray* create_handle(T* t) {
  MexAutoPtr<T>* handle= new MexAutoPtr<T>(t);
  return handle->to_mex_handle();
}

// Obtain object represented by handle.
template <typename T> T& get_object(const mxArray* mxh) {
  MexAutoPtr<T>* handle= MexAutoPtr<T>::from_mex_handle(mxh);
  return handle->get_object();
}

// Obtain indexed object by handle (array of object handles)
template <typename T> T& get_multi_object(const mxArray* mxh, size_t index) {
  
  if (mxGetClassID(mxh) != mxPOINTER_CLASS
      || mxIsComplex(mxh) || mxGetM(mxh)!=1 || mxGetN(mxh) < (index+1) )
    mexErrMsgTxt("Parameter is not an MexAutoPtr array.");

  if (!(index < mxGetDimensions(mxh)[1]))
    mexErrMsgTxt("Index out bounds on MexAutoPtr array.");

  // We *assume* we can store MexAutoPtr<T> pointer in the mxUINT64 of handle
  MexAutoPtr<T>* obj = *reinterpret_cast<MexAutoPtr<T>**>(mxGetPr(mxh)+index);
  if (!obj) {
    // Check to see we don"t have an invalid pointer
    mexErrMsgTxt("Parameter is NULL. It does not represent" \
		 "an MexAutoPtr object.");
  }
  if (obj->signature != obj) {
    // Check memory has correct signature
    mexErrMsgTxt("Parameter does not represent an MexAutoPtr object.");
  }
#ifdef MEX_TYPE_ID_NAME
  if (strcmp(obj->type, typeid(T).name())!=0) {
    mexPrintf("Given: <%s>, Required: <%s>.\n", obj->type, typeid(T).name());
#else
  if (*(obj->type) != typeid(T)) {
    mexPrintf("Given: <%s>, Required: <%s>.\n", obj->type->name(),
	      typeid(T).name());
#endif
    mexErrMsgTxt("Given MexAutoPtr does not represent the correct type.");
  }
  return obj->get_object();
}


// If deleting object, rather than leaving it to garbage collection,
// must delete it via the handle; do not delete T* directly.
template <typename T> void destroy_object(const mxArray *mxh) {
  MexAutoPtr<T>* handle= MexAutoPtr<T>::from_mex_handle(mxh);
  delete handle;
}

// Object generator used for constructing STL containers of objects
// referenced by MexAutoPtr's
template<typename T> class object_generator {
public:
  object_generator(const mxArray* objects) : m_objects(objects), index(0) {
  }
  T operator()() const {
    T object = get_multi_object<T>(m_objects,index);
    index++;
    return object;
  }
  private:
  const mxArray* m_objects;
  mutable size_t index;
};

/*///////////////////////////////////////////////////////////////////
//
// Name:
//   GarbageCollector
//
// Decription:
//
// - Garbage collection singleton (one collector object for each type
//   T). Ensures that registered handles are deleted when the dll is
//   released (they may also be deleted previously without problem).
//   
// - The GarbageCollector provides protection against resource leaks
//   in the case where 'clear all' is called in MatLab. (This is
//   because MatLab will call the destructors of statically allocated
//   objects but not heap-allocated objects.)  Object wrapper
//   class MexAutoPtr is similar to STL auto_ptr
//
///////////////////////////////////////////////////////////////////*/

template <typename T>
class GarbageCollector {
  std::list<MexAutoPtr<T>*> objlist;
  char objname[256];
 public:

  ~GarbageCollector() {
    size_t nObjectsCleared = 0;
    
    typename std::list<MexAutoPtr<T>*>::iterator i;
    typename std::list<MexAutoPtr<T>*>::iterator end = objlist.end();
    for (i= objlist.begin(); i!=end; ++i) {
      // check for valid signature
      if ((*i)->signature == *i) {
        delete *i;
        nObjectsCleared++;
      }
    }
    if (nObjectsCleared) {
#if 0
      mexPrintf("Garbage collector: ");
      mexPrintf("Cleared %zu %s item(s)\n",
		nObjectsCleared, this->objname);
#endif
    }
  }

  static void register_handle (MexAutoPtr<T>* obj) {
    static GarbageCollector singleton(obj);
    singleton.objlist.push_back(obj);
  }

 private: // prevent construction
  GarbageCollector(MexAutoPtr<T>* obj) {
    char buffer[128];
#ifdef MEX_TYPE_ID_NAME
    sprintf(buffer, "%zu", strlen(obj->type));
#else
    sprintf(buffer, "%zu", strlen(obj->type.name()));
#endif

#if defined(__GNUC__)
    strcpy(objname,obj->type+strlen(buffer));
#elif defined(_MSC_VER)
    strcpy(objname,obj->type+6);
#endif

  }
  GarbageCollector(const GarbageCollector&);
};

// Import a handle from MatLab as a mxArray of UINT64. Check that
// it is actually a pointer to an MexAutoPtr<T>.
template <typename T>
MexAutoPtr<T>* MexAutoPtr<T>::from_mex_handle(const mxArray* handle)
{
  if (mxGetClassID(handle) != mxPOINTER_CLASS
      || mxIsComplex(handle) || mxGetM(handle)!=1 || mxGetN(handle)!=1)
    mexErrMsgTxt("Parameter is not an MexAutoPtr type.");

  // We *assume* we can store MexAutoPtr<T> pointer in the mxUINT64 of handle
  MexAutoPtr* obj = *reinterpret_cast<MexAutoPtr<T>**>(mxGetPr(handle));

  if (!obj) {
    // Check to see we don"t have an invalid pointer
    mexErrMsgTxt("Parameter is NULL. It does not represent an" \
		 "MexAutoPtr object.");
  }

  // TODO: change this for max-min check for pointer values

  if (obj->signature != obj) {
    // Check memory has correct signature
    mexErrMsgTxt("Parameter does not represent an MexAutoPtr object.");
  }
#ifdef MEX_TYPE_ID_NAME
  if (strcmp(obj->type, typeid(T).name())!=0) {
    mexPrintf("Givem: <%s>, Required: <%s>.\n", obj->type, typeid(T).name());
#else
  if (*(obj->type) != typeid(T)) {
    mexPrintf("Given: <%s>, Required: <%s>.\n", obj->type->name(),
	      typeid(T).name());
#endif
    mexErrMsgTxt("Given MexAutoPtr does not represent the correct type.");
  }

  return obj;
}

// Create a numeric array as handle for an MexAutoPtr.
// We ASSUME we can store object pointer in the mxUINT64 element of mxArray.
template <typename T>
mxArray* MexAutoPtr<T>::to_mex_handle()  {

  mxArray* handle  = mxCreateNumericMatrix(1, 1, mxPOINTER_CLASS , mxREAL);

  *reinterpret_cast<MexAutoPtr<T>**>(mxGetPr(handle)) = this;
  return handle;
}

#endif
