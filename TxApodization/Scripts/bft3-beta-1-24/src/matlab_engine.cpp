#include "matlab_engine.h"

cMatlabEngine::cMatlabEngine() {
  pEng = NULL;
}
cMatlabEngine::~cMatlabEngine() {
  if (pEng)
    close();
}

int cMatlabEngine::open(const char* cmd) {
  pEng = engOpen(cmd);
  return 0;
}

int cMatlabEngine::close() {
  int result = engClose(pEng);
  if (result==0)
    pEng = NULL;

  return result;
}

int cMatlabEngine::eval_string(const char* str) {
  return (engEvalString(pEng, str));
}

mxArray* cMatlabEngine::get_variable(const char *name){
  return (engGetVariable(pEng, name));
}
int cMatlabEngine::set_variable(const char *name, const mxArray *mp) {
  return (engPutVariable(pEng, name, mp));
}

int cMatlabEngine::register_buffer(char *p, int n) {
  return (engOutputBuffer(pEng, p, n));
}

