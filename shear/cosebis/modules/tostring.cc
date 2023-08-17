#include "tostring.h"

void noBlanks(string& str)
{
 unsigned n;
 for(n=0;n<str.size();n++)
  if ((str[n]==32)||(str[n]==10)) str.erase(n--,1);
}

string toString(const int i){
  ostrstream ost;
  ost << i << '\0';
  return ost.str();
}

string toString(const long l){
  ostrstream ost;
  ost << l << '\0';
  return ost.str();
}

string toString(const int short us){
  ostrstream ost;
  ost << us << '\0';
  return ost.str();
}


string toString(const int unsigned ui){
  ostrstream ost;
  ost << ui << '\0';
  return ost.str();
}

string toString(const float f,int precision){
  ostrstream ost;
 
  if (precision)
    ost.precision(precision);
  else if (precision==-1)
    ost.precision(0);

  ost.setf(ios::fixed,ios::floatfield);
  ost << f << '\0';
  return ost.str();
}

string toString(const double d,int precision){
  ostrstream ost;

  if (precision)
    ost.precision(precision);
  else if (precision==-1)
    ost.precision(0);

  ost.setf(ios::fixed,ios::floatfield);
  ost << d << '\0';
  return ost.str();
}

