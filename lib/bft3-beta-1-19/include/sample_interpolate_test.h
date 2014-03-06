#ifndef SAMPLE_INTERPOLATE_TEST_H
#define SAMPLE_INTERPOLATE_TEST_H

#include <cppunit/extensions/HelperMacros.h>

#include <iostream>
using std::cout;
using std::endl;

#include <stdexcept>
using std::runtime_error;

#include "matlab_engine.h"

class SampleInterpolateTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( SampleInterpolateTest );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testNearestNeighbour );
  CPPUNIT_TEST( testLinear );
  CPPUNIT_TEST( testCubic );
  CPPUNIT_TEST( testSpline );
  CPPUNIT_TEST_EXCEPTION( testEvalThrow, std::runtime_error );
  CPPUNIT_TEST_EXCEPTION( testConstructorThrow1, std::runtime_error );
  CPPUNIT_TEST_SUITE_END();
  
 public:
  void setUp();
  void tearDown();

  void testConstructor();
  void testNearestNeighbour();
  void testLinear();
  void testCubic();
  void testSpline();
  void testConstructorThrow1();
  void testConstructorThrow2();
  void testConstructorThrow3();
  void testConstructorThrow4();
  void testConstructorThrow5();
  void testEvalThrow();

 private:

  Engine* m_eng;
};

#endif 
