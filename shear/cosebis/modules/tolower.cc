#include <cstring>
#include "tolower.h"

void toLower(string& str){
  int s=str.size(),i=0;
  for (i=0;i<s;i++)
    str[i]=tolower(str[i]);
}
