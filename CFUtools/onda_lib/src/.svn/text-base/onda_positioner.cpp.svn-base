/**
 * onda_positioner.cpp
 * Defines the functions that control the positioner of the Onda system.
 * All commands are send as text strings via TCP.
 *
 * $Id: onda_positioner.cpp 6 2012-04-16 15:40:01Z mf $
 */

#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <string>

#include <mex.h>
#include "onda_positioner.h"
#include <tcp.h>
#include "err_printf.h"

using std::string;
using std::ostringstream;
using std::stringstream;


int SetPositionerLowLimit(int axis, float value)
{
  ostringstream ostm;
  ostm << "Positioner " << "Axis" << axis << " LowLimit " << value <<"\n";
  return tcp_tx(ostm.str().c_str());
}


int GetPositionerLowLimit(int axis, float* result)
{
  ostringstream ostm;
  stringstream iostm (stringstream::in | stringstream::out);
  int retval;
  char* rx_buf;

  /* make string */
  ostm << "Positioner " << "Axis" << axis << " LowLimit?\n";
  /* send string */
  retval = tcp_query(ostm.str().c_str(), &rx_buf);
  if (retval){
    err_printf ("Error in GetPositionerLowLimit.\n");
    return -1;
  }

  printf("RX: %s\n", rx_buf);
  /* string to float */
  iostm << rx_buf;
  iostm >> *result;

  return 0;
}

int SetPositionerHighLimit(int axis, float  value)
{
  ostringstream ostm;
  ostm << "Positioner " << "Axis" << axis << " HighLimit " << value <<"\n";
  return tcp_tx(ostm.str().c_str());
}

int GetPositionerHighLimit(int axis, float* result)
{
  ostringstream ostm;
  stringstream iostm (stringstream::in | stringstream::out);
  int retval;
  char* rx_buf;
  char  str_out[50];
  /* make string */
  sprintf(str_out, "Positioner Axis %i HighLimit?\n", axis);
  //ostm << "Positioner " << "Axis" << axis << " HighLimit?\n";
  /* send and receive */
  //retval = tcp_query(ostm.str().c_str(), &rx_buf);
  retval = tcp_query(str_out, &rx_buf);
  if (retval){
    err_printf ("Error in GetPositionerHighLimit.\n");
    return -1;
  }

  printf("RX: %s\n", rx_buf);
  
  /* string to int */
  
  sscanf(rx_buf, "%f",result);
  //iostm << *rx_buf;
  //iostm >> *result;
  printf("after conversion: %f\n", *result);

  return 0;
}

int GetPositionerStepsPerSecond(int axis, int* result)
{
  ostringstream ostm;
  stringstream iostm (stringstream::in | stringstream::out);
  int retval;
  char* rx_buf;
  /* make string */
  ostm << "Positioner " << "Axis" << axis << " StepsPerSecond?\n";

  /* send and receive */
  retval = tcp_query(ostm.str().c_str(), &rx_buf);
  if (retval){
    err_printf ("Error in GetPositionerLowLimit.\n");
    return -1;
  }

  /* string to float */
  iostm << *rx_buf;
  iostm >> *result;

  return 0;
}


int GetPositionerMinStepsPerSecond(int axis, int* result)
{
  ostringstream ostm;
  stringstream iostm (stringstream::in | stringstream::out);
  int retval;
  char* rx_buf;

  /* make string */
  ostm << "Positioner " << "Axis" << axis << " MinStepsPerSecond?";

  /* send and receive */
  retval = tcp_query(ostm.str().c_str(), &rx_buf);
  if (retval){
    err_printf ("Error in GetPositionerMinStepsPerSecond.\n");
    return -1;
  }

  /* string to float */
  iostm << *rx_buf;
  iostm >> *result;

  return 0;
}



int  PositionerMoveRel(int axis, float  value)
{
  ostringstream ostm;
  ostm << "Positioner " << "Axis" << axis << " MoveRel " << value <<"\n";
  return tcp_tx(ostm.str().c_str());
}

int PositionerMoveAbs(int axis, float value)
{
  ostringstream ostm;
  ostm << "Positioner " << "Axis" << axis << " MoveAbs " << value <<"\n";
  return tcp_tx(ostm.str().c_str());
}

int SetPosition(int axis, float value)
{
  ostringstream ostm;
  ostm << "Positioner " << "Axis" << axis << " Position " << value <<"\n";
  return tcp_tx(ostm.str().c_str());
}



int GetPosition(int axis,float* result)
{
  ostringstream ostm;
  stringstream iostm (stringstream::in | stringstream::out);
  int retval;
  char* rx_buf;

  *result= 42421;
  printf("In GetPosition1\n");
  /* make string */
  ostm << "Positioner " << "Axis" << axis << " Position?\n";
  
  info_printf("In GetPosition2\n");
  info_printf("sending: %s\n", ostm.str().c_str());

  /* send and receive */
  retval = tcp_query(ostm.str().c_str(), &rx_buf);
  info_printf("retval: %i\n", retval);
  if (retval){
      err_printf ("Error sending in GetPosition.\n");
    return -1;
  }
  info_printf("received: %s.\n", rx_buf);
  info_printf("In GetPosition3\n");
  
  
  //sscanf(rx_buf, "%f", result);
  //info_printf("converted: %f\n", *result);
  
  /* string to float */
  //iostm << rx_buf;
  //iostm >> *result;
  printf("In GetPosition4\n");

  return 0;
}










/*

int SetPositionerStepsPerSecond(axis, int value		      )
{
  tcp_tx("Positioner " << "Axis" << axis << " StepsPerSecond " << IntToStr(value));
}

int SetPositionerMinStepsPerSecond(axis, int value			 )
{
  tcp_tx("Positioner " << "Axis" << axis << " MinStepsPerSecond " << IntToStr(value));
}

int GetPositionerAxisCount: integer; export; stdcall;
{
  result = QueryInt("Positioner AxisCount?");
}

int GetPositionerReversed(int axis): longbool; export; stdcall;
{
  result = QueryBool("Positioner " << "Axis" << axis << " Reversed?");
}

int SetPositionerReversed(int axis, value: longbool); export; stdcall;
{
  tcp_tx("Positioner " << "Axis" << axis << " Reversed " << BoolToStr(value));
}

int GetUseMachineCoordinates: longbool; export; stdcall;
{
  result = QueryBool("Positioner UseMachineCoordinates?");
}

int SetUseMachineCoordinates(value: longbool); export; stdcall;
{
  tcp_tx("Positioner UseMachineCoordinates " << BoolToStr(value));
}








int FindLimitSwitch(int axis, value: longbool)
var
  s: string;
{
  if value then s = "Positive" else s = "Negative";
  tcp_tx("Positioner " << "Axis" << axis << " FindLimitSwitch " << s);
}

int GetLimitSwitchState(int axis, value: longbool	     )
{
  result = QueryBool("Positioner " << "Axis" << axis << " LimitSwitch? " << BoolToStr(value));
}






int GetPositionerRampSteps(int axis		)
{
  result = QueryInt("Positioner " << "Axis" << axis << " RampSteps?");
}


int SetPositionerRampSteps(axis, int value		 )
{
  tcp_tx("Positioner " << "Axis" << axis << " RampSteps " << IntToStr(value));
}



int GetPositionerModel(char* result	     )
{
  QueryStr("Positioner Model?", value);
}

int GetPositionerSerial(char* result	      )
{
  QueryStr("Positioner Serial?", value);
}

int GetPositionerFirmwareVersion(char* result		)
{
  QueryStr("Positioner FirmwareVersion?", value);
}

int GetPositionerTranslation(int axis)
{
  result = QueryBool("Positioner " << "Axis" << axis << " Translation?");
}

int SetPositionerTranslation(int axis, value: longbool); export; stdcall;
{
  tcp_tx("Positioner " << "Axis" << axis << " Translation " << BoolToStr(value));
}

int GetPositionerAxisName(int axis, char* result); export; stdcall;
{
  QueryStr("Positioner " << "Axis" << axis << " Name?", value);
}

int GetPositionerUnits(int axis, char* result); export; stdcall;
{
  QueryStr("Positioner " << "Axis" << axis << " Units?", value);
}

int GetPositionerScaleFactor(int axis): double; export; stdcall;
{
  result = QueryDouble("Positioner " << "Axis" << axis << " ScaleFactor?");
}


int GetAngPosInstalled: longbool; export; stdcall;
{
  result = QueryBool("Positioner AngPosInstalled?");
}

int GetAngPosOrientation: integer; export; stdcall;
{
  result = QueryInt("Positioner AngPosOrientation?");
}

int GetAngPosArmOrientation: integer; export; stdcall;
{
  result = QueryInt("Positioner AngPosArmOrientation?");
}

*/
