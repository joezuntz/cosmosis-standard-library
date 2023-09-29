#include "tostring.h"

void noBlanks(string& str)
{
 unsigned n;
 for(n=0;n<str.size();n++)
  if ((str[n]==32)||(str[n]==10)) str.erase(n--,1);
}

string toString(const int i){
  ostringstream ost;
  ost << i << '\0';
  return ost.str();
}

string toString(const long l){
  ostringstream ost;
  ost << l << '\0';
  return ost.str();
}

string toString(const int short us){
  ostringstream ost;
  ost << us << '\0';
  return ost.str();
}


string toString(const int unsigned ui){
  ostringstream ost;
  ost << ui << '\0';
  return ost.str();
}

string toString(const float f,int precision){
  ostringstream ost;
 
  if (precision)
    ost.precision(precision);
  else if (precision==-1)
    ost.precision(0);

  ost.setf(ios::fixed,ios::floatfield);
  ost << f << '\0';
  return ost.str();
}

string toString(const double d,int precision){
  ostringstream ost;

  if (precision)
    ost.precision(precision);
  else if (precision==-1)
    ost.precision(0);

  ost.setf(ios::fixed,ios::floatfield);
  ost << d << '\0';
  return ost.str();
}

