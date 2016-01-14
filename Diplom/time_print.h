// TIME_PRINT Print out the current date and time.
// Notes:
// The standard Fortran 90 routine DATE_AND_TIME is used to get
// the current date and time strings.

#include <ctime>
#include <iostream>

using namespace std;

//void Time_Print()
//{
//	char buffer[80];
//	time_t seconds = time(NULL);
//	tm* timeinfo = localtime(&seconds);
//	char* format = "%A, %B %d, %Y %I:%M:%S";
//	strftime(buffer, 80, format, timeinfo);
//	cout << "Current Datetime: " << buffer << endl;
//
//
//  
//  //const int output = 6;
//  ////local scalars
//  //character ( len = 8 ) :: datstr
//  //character ( len = 10 ) :: timstr
//  //// Get the current date and time
//  //date_and_time ( datstr, timstr )
//  //// Write out the date and time.
//  //write ( output, "(/A)" ) " Date = " // datstr(7:8) //"/" // &
//  //	datstr(5:6) // "/" // &
//  //	datstr(1:4)
//  //write ( output, "(A)" ) " Time = " // timstr(1:2) //":" // &
//  //	timstr(3:4) // ":" // &
//  //	timstr(5:10)
//  //write ( output, *)
//}
