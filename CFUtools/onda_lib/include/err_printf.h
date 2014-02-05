
#ifndef C_PRINT_H
#define C_PRINT_H

#ifndef MATLAB_MEX_FILE
#define err_printf   printf
#define info_printf  printf
#else
#include "mex.h"
#define info_printf  mexPrintf
#define err_printf   mexErrMsgTxt
#endif

#endif

