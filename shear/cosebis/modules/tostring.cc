#include "tostring.h"
#include <string>
#include <sstream>

void noBlanks(string& str)
{
 unsigned n;
 for(n=0;n<str.size();n++)
  if ((str[n]==32)||(str[n]==10)) str.erase(n--,1);
}

string toString(const int i){
  return to_string(i);
}

string toString(const long l){
  return to_string(l);
}

string toString(const int short us){
  return to_string(us);
}


string toString(const int unsigned ui){
  return to_string(ui);
}

string toString(const float f,int precision){
  ostringstream ost;
 
  if (precision)
    ost.precision(precision);
  else if (precision==-1)
    ost.precision(0);

  ost.setf(ios::fixed,ios::floatfield);
  ost << f;
  return ost.str();
}

string toString(const double d,int precision){
  ostringstream ost;

  if (precision)
    ost.precision(precision);
  else if (precision==-1)
    ost.precision(0);

  ost.setf(ios::fixed,ios::floatfield);
  ost << d;
  return ost.str();
}

