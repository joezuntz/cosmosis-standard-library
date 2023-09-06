#include "ticker.h"

void ticker::reset()
{
  time_t Time=time(&Time);
  started = (int) Time;
}

ticker::ticker()
{
  reset();
  clog.precision(2);
}

ticker::~ticker()
{
  clog << "\r\r";
  clog.flush();
}

void ticker::tick(number percentage)
{
  static float diff;
  static char  unit;

  time_t Time=time(&Time);
  diff = (1.-percentage)/percentage*(Time-started);

  if (diff<minute)
    unit=0;
  else if (diff<hour)
    {unit=1;diff/=minute;}
  else if (diff<day)
    {unit=2;diff/=hour;}
  else if (diff<week)
    {unit=3;diff/=day;}
  else if (diff<year)
  {unit=4;diff/=week;}
  else if (diff<century)
    {unit=5;diff/=year;}
  else
    {unit=6;diff/=century;}

  clog << "Progress: " 
       << (100.0*percentage);

  switch (unit)
   {
   case 0 : clog << "sec]";break;
   case 1 : clog << "min]";break;
   case 2 : clog << "h]";break;
   case 3 : clog << "days]";break;
   case 4 : clog << "weeks]";break;
   case 5 : clog << "yrs]";break;
   case 6 : clog << "cent]";break;
   }

  clog << "                        \r\r";
  clog.flush();
}

number ticker::runtime()
{
  time_t Time=time(&Time);
  int diff = (int) Time-started;

  return diff/60.0;
}
