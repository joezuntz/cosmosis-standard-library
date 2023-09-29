#ifndef TOSTRING_H
#define TOSTRING_H

#include <string>

using namespace std;

/** Original code by Patrick Simon
 * 
 * Trunctuates all blanks in front and in back of string 
 */
void noBlanks(string& str); 	

/** Converts long to string
*/
string toString(const long);

/** Converts int to string
 */
string toString(const int);

/** Converts short to string
 */
string toString(const int short);

/** Converts int to string
 */
string toString(const unsigned);

/** Converts foat to string
 */
string toString(const float,int=0);

/** Converts double to string
 */
string toString(const double,int=0);

#endif
