#ifndef APODIZATION_TEST_H
#define APODIZATION_TEST_H

#include <cppunit/extensions/HelperMacros.h>

class ApodizationTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( ApodizationTest );
  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST_SUITE_END();
  
 public:
  void setUp();
  void tearDown();

  void testConstructor();
 private:

};

#endif 
