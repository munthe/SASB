// $Id: sample_interpolate_test.cpp,v 1.12 2011-05-01 21:25:56 jmh Exp $
#include "sample_interpolate_test.h"
#include "sample_interpolate.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION( SampleInterpolateTest );

void 
SampleInterpolateTest::setUp() {
}

void 
SampleInterpolateTest::tearDown() {
}

void 
SampleInterpolateTest::testConstructor() {
}

void 
SampleInterpolateTest::testNearestNeighbour() {
  cMatlabEngine matlab; 

  // Open engine
  if (matlab.open("\0") != 0) {
    return;
  }
  
  const char* s1 =
    "fs=20;"
    "f0=5;"
    "s_offset = 2.3;"
    ""
    "rand('seed',0.23);"
    ""
    "% Signal\n"
    "s = (0:1:2*fs/f0);"
    "v = sin(2*pi*f0/fs*(s+s_offset));"
    ""
    "% Random sample locations\n"
    "random_samples = (2*fs/f0)*[rand(1,5)];"
    ""
    "% Nearest neighbour\n"
    "nn_v = interp1(s,v,random_samples,'nearest');";
  
  if (matlab.eval_string(s1)) {
    cout << "Error executing Matlab expression" << endl;
    return;
  }

  mxArray* mxValues = matlab.get_variable("v");
  mxArray* mxRandomSamples = matlab.get_variable("random_samples");
  mxArray* mxNearestNeighbor = matlab.get_variable("nn_v");
  
  mwSize n_samples, n_random_samples;

  // Samples
  n_samples = mxGetDimensions(mxValues)[1];
  const double* samples = (const double*) mxGetData(mxValues);

  // Random locations
  n_random_samples  = mxGetDimensions(mxRandomSamples)[1];
  const double* rnd_locations = (const double*) mxGetData(mxRandomSamples);

  // Nearest neighbor reference
  const double* nearest  =  (const double*) mxGetData(mxNearestNeighbor);

  // Initialize interpolator
  SampleInterpolate<size_t,double> interp = SampleInterpolate<size_t,double>(samples,n_samples);
  interp.setMethod(SampleInterpolate<size_t,double>::nearestNeighbour);

  double difference = 0.0;

  for (size_t i = 0; i < n_random_samples ; i++) {
    difference =std::max(difference,fabs(interp(rnd_locations[i]) - nearest[i]));;
  }

  mxDestroyArray(mxValues);
  mxDestroyArray(mxRandomSamples);
  mxDestroyArray(mxNearestNeighbor);

  CPPUNIT_ASSERT( difference < DBL_EPSILON);
}

void 
SampleInterpolateTest::testLinear() {
  cMatlabEngine matlab; 

  // Open engine
  if (matlab.open("\0") != 0) {
    return;
  }
  
  const char* s1 =
    "fs=20;"
    "f0=5;"
    "s_offset = 2.3;"
    ""
    "rand('seed',0.23);"
    ""
    "% Signal\n"
    "s = (0:1:2*fs/f0);"
    "v = sin(2*pi*f0/fs*(s+s_offset));"
    ""
    "% Random sample locations\n"
    "random_samples = (2*fs/f0)*[rand(1,5)];"
    ""
    "% Linear\n"
    "l_v = interp1(s,v,random_samples,'linear');";

  if (matlab.eval_string(s1)) {
    cout << "Error executing Matlab expression" << endl;
    return;
  }

  mxArray* mxValues = matlab.get_variable("v");
  mxArray* mxRandomSamples = matlab.get_variable("random_samples");
  mxArray* mxLinear = matlab.get_variable("l_v");
  
  mwSize n_samples, n_random_samples;

  // Samples
  n_samples = mxGetDimensions(mxValues)[1];
  const double* samples = (const double*) mxGetData(mxValues);

  // Random locations
  n_random_samples  = mxGetDimensions(mxRandomSamples)[1];
  const double* rnd_locations = (const double*) mxGetData(mxRandomSamples);

  // Linear interpolation reference
  const double* linear  =  (const double*) mxGetData(mxLinear);

  // Initialize interpolator
  SampleInterpolate<size_t,double> interp = SampleInterpolate<size_t,double>(samples,n_samples);

  double difference = 0.0;

  for (size_t i = 0; i < n_random_samples ; i++) {
    difference =std::max(difference,fabs(interp(rnd_locations[i]) - linear[i]));;
  }

  CPPUNIT_ASSERT( difference < DBL_EPSILON);

  mxDestroyArray(mxValues);
  mxDestroyArray(mxRandomSamples);
  mxDestroyArray(mxLinear);
}

void 
SampleInterpolateTest::testCubic() {
  cMatlabEngine matlab; 

  // Open engine
  if (matlab.open("\0") != 0) {
    return;
  }
  
  const char* s1 =
    "fs = 20;"   
    "f0=5;"
    "s_offset = 2.3;"
    ""
    "rand('seed',0.23);"
    ""
    "% Signal\n"
    "s = (0:1:2*fs/f0);"
    "v = sin(2*pi*f0/fs*(s+s_offset));"
    ""
    "% Random sample locations\n"
    "random_samples = (2*fs/f0)*[rand(1,5)];"
    ""
    "% Polynomial\n"
    "p = [];"
    "p_t = [];"
    ""
    "p_s = [];"
    "p_t_s = [];"
    ""
    "pp = polyfit(s(1:4), v(1:4),3);"
    ""
    "inx = find(random_samples < s(3));"
    "if ~isempty(inx)"
    "  ss = random_samples(inx);"
    "  val = pp(1)*ss.^3 + pp(2)*ss.^2 + pp(3)*ss + pp(4);"
    "  p_t_s = [p_t_s ss];"
    "  p_s = [p_s val];"
    "end;"
    ""
    "for i=3:length(v)-3"
    "  pp = polyfit(s(i:i+3),v(i:i+3),3);"
    "  inx = find((random_samples>=s(i)) & (random_samples < s(i+1)));"
    "  if ~isempty(inx)"
    "    ss = random_samples(inx);"
    "    p_t_s = [p_t_s ss];"
    "    val = pp(1)*ss.^3 + pp(2)*ss.^2 + pp(3)*ss + pp(4);"
    "    p_s = [p_s val];"
    "  end;"
    "end;"
    ""
    "inx = find(random_samples >= s(length(v)-2));"
    "pp = polyfit(s(length(v)-3:length(v)), v(length(v)-3:length(v)),3);"
    "if ~isempty(inx)"
    "  ss = random_samples(inx);"
    "  val = pp(1)*ss.^3 + pp(2)*ss.^2 + pp(3)*ss + pp(4);"
    "  p_t_s = [p_t_s ss];"
    "  p_s = [p_s val];"
    "end;";

  if (matlab.eval_string(s1)) {
    cout << "Error executing Matlab expression" << endl;
    return;
  }
  
  mxArray* mxValues = matlab.get_variable("v");
  mxArray* mxRandomSamples = matlab.get_variable("random_samples");
  mxArray* mxCubic = matlab.get_variable("p_t_s");
  
  mwSize n_samples, n_random_samples;

  // Samples
  n_samples = mxGetDimensions(mxValues)[1];
  const double* samples = (const double*) mxGetData(mxValues);

  // Random locations
  n_random_samples  = mxGetDimensions(mxRandomSamples)[1];
  const double* rnd_locations = (const double*) mxGetData(mxRandomSamples);

  const double* cubic = (const double*) mxGetData(mxCubic);

  // Initialize interpolator
  SampleInterpolate<size_t,double> interp = SampleInterpolate<size_t,double>(samples,n_samples);
  interp.setMethod(SampleInterpolate<size_t,double>::cubic);

  double difference = 0.0;

  for (size_t i = 0; i < n_random_samples ; i++) {
    difference =std::max(difference,fabs(interp(rnd_locations[i]) - cubic[i]));
    cout << "Difference: " << difference << endl;
  }

  mxDestroyArray(mxValues);
  mxDestroyArray(mxRandomSamples);
  mxDestroyArray(mxCubic);

  cout << "Difference: " << difference << endl;
  CPPUNIT_ASSERT( difference < DBL_EPSILON);

}

void 
SampleInterpolateTest::testSpline() {
  cMatlabEngine matlab; 

  // Open engine
  if (matlab.open("\0") != 0) {
    return;
  }
  
  const char* s1 =
    "fs=20;"
    "f0=5;"
    "s_offset = 2.3;"
    ""
    "rand('seed',0.23);"
    ""
    "% Signal\n"
    "s = (0:1:2*fs/f0);"
    "v = sin(2*pi*f0/fs*(s+s_offset));"
    ""
    "% Random sample locations\n"
    "random_samples = (2*fs/f0)*[rand(1,5)];"
    ""
    "% Spline fit 2nd derivative continuous\n"
    "cs_v = interp1(s,v,random_samples,'spline');";

  if (matlab.eval_string(s1)) {
    cout << "Error executing Matlab expression" << endl;
    return;
  }

  mxArray* mxValues        = matlab.get_variable("v");
  mxArray* mxRandomSamples = matlab.get_variable("random_samples");
  mxArray* mxCubicSpline   = matlab.get_variable("cs_v");
  
  mwSize n_samples, n_random_samples;

  // Samples
  n_samples = mxGetDimensions(mxValues)[1];
  const double* samples = (const double*) mxGetData(mxValues);

  // Random locations
  n_random_samples  = mxGetDimensions(mxRandomSamples)[1];
  const double* rnd_locations = (const double*) mxGetData(mxRandomSamples);

  // Linear interpolation reference
  const double* cubic_spline  =  (const double*) mxGetData(mxCubicSpline);

  // Initialize interpolator
  SampleInterpolate<size_t,double> interp = SampleInterpolate<size_t,double>(samples,n_samples);
  interp.setMethod(SampleInterpolate<size_t,double>::spline);

  double difference = 0.0;
  double interpolated_value;

  for (size_t i = 0; i < n_random_samples ; i++) {
    interpolated_value = interp(rnd_locations[i]);
    difference =std::max(difference,fabs(interpolated_value - cubic_spline[i]));
  }

  mxDestroyArray(mxValues);
  mxDestroyArray(mxRandomSamples);
  mxDestroyArray(mxCubicSpline);

//  CPPUNIT_ASSERT( difference < DBL_EPSILON);

}

void
SampleInterpolateTest::testEvalThrow() {
  double samples[3] = {0.0,1.0,2.0};
  size_t n_samples = 3;
  SampleInterpolate<size_t,double> interp = SampleInterpolate<size_t,double>(samples,n_samples);
  interp.setMethod(SampleInterpolate<size_t,double>::cubic);
  interp(0.8);
}

void
SampleInterpolateTest::testConstructorThrow1() {
  double* samples  = NULL;
  size_t n_samples = 0;
  SampleInterpolate<size_t,double> interp = SampleInterpolate<size_t,double>(samples,n_samples);
}
