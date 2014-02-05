/**
 * onda_lib.cpp
 * This file makes the connection between MATLAB and the C++ library which
 * connects to the Onda system.
 *
 */




#include <cstdio>
#include <cstring>
#include <mex.h>
#include "onda_positioner.h"
#include "tcp.h"
#include "err_printf.h"

#define HOST_PORT 49999

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  int axis = 0, retval;
  float value = 0;
  float* value_p = NULL;

  /* string buffers */
  enum {host_addr_len  = 20};
  enum {func_str_len   = 30};
  char host_addr[host_addr_len];
  char func_str[func_str_len];

  //  static float atest = 10.1;
  
  /* output pointers */
  mxArray* val_p = NULL;

  if (nrhs < 1) err_printf("Minimum one input required.");
  if (nlhs > 1) err_printf("Maximum one output argument allowed.");

  /* get input vars */
  mxGetString(prhs[0], func_str, func_str_len);
  if (nrhs > 1)  axis  = (int)mxGetScalar(prhs[1]);
  if (nrhs > 2)  value = (float)mxGetScalar(prhs[2]);

  /* set output pointer */
  val_p = plhs[0];
  val_p = mxCreateDoubleMatrix(1, 1, mxREAL);




  /* Get Version of this library */
  // TODO: Get define the version during compiling.
  if (strcmp("version", func_str) == 0){
      info_printf("Onda connection library, version 1.3, 2013-04-22.\n");
  }


  /* move relative */
  else if (strcmp("move_relative", func_str) == 0){
    if (nrhs != 3) err_printf("This function requires three input arguments.");

    retval = PositionerMoveRel(axis, value);
    if (retval){
      err_printf("Error in Positioner move relative\n");
      return;
    }
  }


  /* move absolute */
  else if (strcmp("move_absolute", func_str) == 0){
    if (nrhs != 3) err_printf("This function requires three input arguments.");

    retval = PositionerMoveAbs(axis, value);
    if (retval) err_printf("Error in positioner move absolut.");
  }


  /* set position */
  else if (strcmp("set_position", func_str) == 0){
    if (nrhs != 3) err_printf("This function requires three input arguments.");

    retval = SetPosition(axis, value);
    if (retval) err_printf("Error in Set position.");
  }


  /* get position */
  else if (strcmp("get_position", func_str) == 0){
    info_printf("debugging\n");
    info_printf("debugging2\n");
    retval = GetPosition(axis, value_p);
    info_printf("debugging11\n");
    mxGetPr(val_p)[0] = value;
    if (0) {
        if (nrhs != 2) err_printf("This function requires two input arguments.");
    
        info_printf("get pos\n");
    
        retval = GetPosition(axis, value_p);
        if (retval) err_printf("Error in Get position.");
        info_printf("Onda pos: %f\n", *value_p);
        /* set output */
        mxGetPr(val_p)[0] = *value_p;}
  }


  /* Get low limit */
  else if (strcmp("get_low_limit", func_str) == 0){
    if (nrhs != 2) err_printf("This function requires two input arguments.");

    retval = GetPositionerLowLimit(axis, value_p);
    if (retval) err_printf("Onda: Error in Get low limit.");
    /* set output */
    mxGetPr(val_p)[0] = *value_p;
  }


  /* set low limit */
  else if (strcmp("set_low_limit", func_str) == 0){
    if (nrhs != 3) err_printf("This function requires three input arguments.");
    
    
    retval = SetPositionerLowLimit(axis, value);
    if (retval) err_printf("Onda: Error in Set low limit.");
  }


  /* Get high limit */
  else if (strcmp("get_high_limit", func_str) == 0){
    if (nrhs != 2) err_printf("This function requires two input arguments.");

    retval = GetPositionerHighLimit(axis, value_p);
    if (retval) err_printf("Error in Get high limit.");
    info_printf("Onda: Pos low limit: %f\n", *value_p);
    /* set output */
    mxGetPr(val_p)[0] = *value_p;
  }


  /* Set high limit */
  else if (strcmp("set_high_limit", func_str) == 0){
    if (nrhs != 3) err_printf("This function requires three input arguments.");

    retval = SetPositionerHighLimit(axis, value);
    if (retval) err_printf("Error in Set high limit.");
  }


  /* Init connection */
  else if (strcmp("init_connection", func_str) == 0){
    if (nrhs != 2) err_printf("This function requires two input arguments.");

    mxGetString(prhs[1], host_addr, host_addr_len);
    retval = tcp_init (host_addr, HOST_PORT);
    if (retval)
      err_printf("Onda: Could not initialise connection.\n");
    else
      info_printf("Onda: connection initialised.\n");

    /* automatic close connection when MATLAB clears or exits */
    mexAtExit(tcp_auto_term);
  }


  /* terminate connection */
  else if (strcmp("terminate_connection", func_str) == 0){
    retval = tcp_term ();
    if (retval)
      err_printf("Onda: Could not terminate connection.\n");
    else
      info_printf("Onda connection terminated.\n");
  }




  /*  CMD not found*/
  else {
    err_printf("Command not found.");
  }
}
