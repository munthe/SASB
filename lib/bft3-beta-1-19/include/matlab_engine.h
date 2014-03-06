#ifndef MATLAB_ENGINE_H
#define MATLAB_ENGINE_H

#include <cstdlib>
#include <cstdio>
#include "engine.h"

class cMatlabEngine {
 public:
  cMatlabEngine();
  ~cMatlabEngine();

  int open(const char* cmd);
  int close();

  int eval_string(const char* str);

  mxArray* get_variable(const char *name);
  int set_variable(const char *name, const mxArray *mp);
  int register_buffer(char *p, int n);

  /*
  int open(const string &cmd);
  int eval_string(const string &str);
  mxArray* get_variable(const string& name);
  int set_variable(const string& name, const mxArray *mp);
  */
 protected:
  Engine* pEng;
  
 private:
};

#endif
