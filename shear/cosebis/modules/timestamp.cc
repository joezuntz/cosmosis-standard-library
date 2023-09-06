#include <cstdlib>
#include <ctime>
#include <iostream>
#include <math.h>

#include "timestamp.h"
#include "tostring.h"

// --------------------------------------------------------------------------

  timestamp::timestamp(): 
    TimestampString()
  {
    setCurrent();
  }

// --------------------------------------------------------------------------

  timestamp::timestamp(const timestamp& stamp){
    *this=stamp;
  }

// --------------------------------------------------------------------------
  
  timestamp::timestamp(int second, int minute, int hour, int day, int month, int year):
    TimestampString()
  {
    setTimeStamp(second, minute, hour, day, month, year);
  }

// --------------------------------------------------------------------------
    
  timestamp::timestamp(const string& timestampstr):
    TimestampString()
  {
    parseString(timestampstr);
  }

// --------------------------------------------------------------------------

  void timestamp::setTimeStamp(int second, int minute, int hour, int day, int month, int year){
    Second=second;  
    Minute=minute;
    Hour=hour;
    Day=day;
    Month=month;
    Year=year;
    Seconds=calculateTMSeconds();  
    setTimeStampString();
  }
// --------------------------------------------------------------------------

  void timestamp::setTimeStamp(const timestamp& stamp){
    setSeconds(stamp.getSeconds());    
  }

// --------------------------------------------------------------------------

  
  void timestamp::setSeconds(const long int seconds){
    struct tm* TM = localtime(&seconds);
    setMembers(*TM);
  }

// --------------------------------------------------------------------------
  
  void timestamp::parseString(const string& timestampstr){
    // To Do: Errorhandling
    TimestampString=timestampstr;
    string YearStr=TimestampString.substr(0,4);
    Year=atoi(YearStr.c_str())-1900;
    string MonthStr=TimestampString.substr(4,2);
    Month=atoi(MonthStr.c_str())-1;
    string DayStr=TimestampString.substr(6,2);
    Day=atoi(DayStr.c_str());    
   	string HourStr=TimestampString.substr(8,2);
   	Hour=atoi(HourStr.c_str());
   	string MinuteStr=TimestampString.substr(10,2);
   	Minute=atoi(MinuteStr.c_str());
   	string SecondStr=TimestampString.substr(12,2);
   	Second=atoi(SecondStr.c_str());
   	Seconds=calculateTMSeconds();
		
  }

// --------------------------------------------------------------------------

  int timestamp::getSeconds() const{
    return Seconds;
  }
// --------------------------------------------------------------------------
  int timestamp::toSeconds(const string& secstr){
    return atoi(secstr.c_str());
  // To Do: improve implemententaion
  }
// --------------------------------------------------------------------------

  string timestamp::getAsString() const
  {
    return TimestampString;
  }
  
// --------------------------------------------------------------------------

  string timestamp::getFormatted() const{
    string ret("");
    ret+=toString(Day);
    ret+=".";
    ret+=toString(Month+1);
    ret+=".";
    ret+=toString(1900+Year);
    ret+="  ";
    ret+=toString(Hour);
    ret+=":";
    ret+=toString(Minute);
    ret+=":";
    ret+=toString(Second);
    return ret;
  }
  
// --------------------------------------------------------------------------
  void timestamp::setCurrent(){
    time_t Time=time(&Time);
    Seconds=(int)Time;
    struct tm* TM=localtime(&Time);
    setMembers(*TM);
  }
// --------------------------------------------------------------------------
  void timestamp::setMembers(const struct tm& TM)
  {
    struct tm _TM = TM;
    Second=TM.tm_sec;
    Minute=TM.tm_min;
    Hour=TM.tm_hour;
    Day=TM.tm_mday;
    Month=TM.tm_mon;
    Year=TM.tm_year;
    DayOfWeek = TM.tm_wday;
    DayOfYear = TM.tm_yday;
    Seconds=(int)mktime(&_TM);  
    setTimeStampString();
  }

  void timestamp::setTimeStampString(){
    TimestampString=toString(1900+Year);
    if (Month<9)
      TimestampString+="0";
    TimestampString+=toString(Month+1);
    if (Day<10)
      TimestampString+="0";
    TimestampString+=toString(Day);
    if (Hour<10)
      TimestampString+="0";
    TimestampString+=toString(Hour);
    if (Minute<10)
      TimestampString+="0";
    TimestampString+=toString(Minute);
    if (Second<10)
      TimestampString+="0";
    TimestampString+=toString(Second);
  }


// --------------------------------------------------------------------------
  int timestamp::calculateTMSeconds(){
    struct tm TM;
    TM.tm_sec=Second;
    TM.tm_min=Minute;
    TM.tm_hour=Hour;
    TM.tm_mday=Day;
    TM.tm_mon=Month;
    TM.tm_year=Year;
    TM.tm_isdst=-1;

    return (int)mktime(&TM);
  }
// --------------------------------------------------------------------------

  bool operator==(const timestamp& t1, const timestamp& t2){
    return (t1.getSeconds() == t2.getSeconds());
  }
  
  bool operator!=(const timestamp& t1, const timestamp& t2){
    return (t1.getSeconds() != t2.getSeconds());
  }
  
  bool operator<(const timestamp& t1, const timestamp& t2){
    return (t1.getSeconds() < t2.getSeconds());
  }
  
  bool operator<=(const timestamp& t1, const timestamp& t2){
    return (t1.getSeconds() <= t2.getSeconds());
  }
  
  bool operator>(const timestamp& t1, const timestamp& t2){
    return (t1.getSeconds() > t2.getSeconds());
  }
  
  bool operator>=(const timestamp& t1, const timestamp& t2){
    return (t1.getSeconds() >= t2.getSeconds());
  }
  
  timestamp operator+(const timestamp& rhs, const int seconds){
    timestamp result(rhs);
    result.setSeconds(result.getSeconds()+seconds);
    return result;
  }  
  
  timestamp& timestamp::operator=(const timestamp& stamp){
    setTimeStamp(stamp);
    return *this;
  }


string  timestamp::getMonthStr(int land)
{
 if (!land)
  switch (Month){
      case 0: return "Jan.";break;
      case 1: return "Feb.";break;
      case 2: return "March";break;
      case 3: return "April"; break;
      case 4: return "May"; break;
      case 5: return "June";break;
      case 6: return "July";break;
      case 7: return "Aug.";break;
      case 8: return "Sep.";break;
      case 9: return "Oct.";break;
      case 10: return "Nov.";break;
      case 11: return "Dec.";break;
      default:
         return "Fehler in der Datumsroutine";
        break;
    }
 else if (land==1)
  switch (Month){
      case 0: return "Januar";break;
      case 1: return "Feburar";break;
      case 2: return "M&auml;rz";break;
      case 3: return "April"; break;
      case 4: return "Mai"; break;
      case 5: return "Juni";break;
      case 6: return "Juli";break;
      case 7: return "August";break;
      case 8: return "September";break;
      case 9: return "Oktober";break;
      case 10: return "November";break;
      case 11: return "Dezember";break;
      default:
         return "Fehler in der Datumsroutine";
        break;
    }

 return "noname";
}

  string timestamp::getNiceForm(int land) const
{
    string ret("");

    if (Day<10) ret+="0";
	
    ret+=toString(Day);
    ret+=".";
    if (Month<9)
      ret+="0";
    ret+=toString(Month+1);// getMonthStr(land);
    ret+=".";

  return (ret+toString(Year+1900));
}

double timestamp::getJulianDate()
{
 return julian(Year+1900,Month+1,Day,Hour,Minute,Second);
}

//---- taken from the book Duffet-Smith, Peter, 
//---- Astronomy with your Personal Computer, Cambridge University Press


double sgn  (double x){
 if (x<0) return -1; 
  else if (x>0) return +1;

 return 0;}
double fnlif(double x){return floor(x);}
double fnunw(double x,double y){return x-fnlif(x/y)*y;}
double fnrad(double x){return 1.745329252E-2*x;}
double fndeg(double x){return 5.729577951E1*x;}

double julian(double year,double month, double day, 
             double hour,double minute,double second)
{
 double A,B,C,D;

 if (month<=2.0) {month+=12.0;year-=1.0;}
 A = floor(year/100.0);
 if ((year<=1582.0)&&(month<=10.0)&&(day<=15)) B = 0.0;
    else B = 2.0-A+floor(A/4.0);
 C = floor(365.25*year);
 D = floor(30.6001*(month+1.0));

 return (B+C+D+day+1720994.5+hour/24.0+minute/1440.0+second/86400.0);
} 

double timestamp::nearestNewMoon()
{
 double  djd0, djd, k, tn, tf, t, t2,b,c,ms,mm,f,ddjd,a1,b1,a;
 
 djd0 = julian(Year+1900,1,0,0,0,0);
 djd  = julian(Year+1900,Month+1,Day,0,0,0);
 k    = fnlif((Year+(djd-djd0)/365.0)*12.3685+0.5);
 tn   = k/1236.85;
 tf   = (k+0.5)/1236.85;
 t    = tn;

 t2 = t*t;a=29.53*k;
 c  = fnrad(166.56+(132.87-9.173E-3*t)*t);
 b  = 5.8868E-4*k+(1.178E-4-1.55E-7*t)*t2
      +3.3E-4*sin(c)+7.5933E-1;
 ms = 359.2242+360.0*fnunw(k/1.236886E1,1.0)
      -(3.33E-5+3.47E-6*t)*t2;
 mm = 306.0253+360.0*fnunw(k/9.330851E-1,1.0)
      +(1.07306E-2+1.236E-5*t)*t2;
 f  = 21.2964+360.0*fnunw(k/9.214926E-1,1.0)
      -(1.6528E-3+2.39E-6*t)*t2;
 ms = fnunw(ms,360.0);
 mm = fnunw(mm,360.0);
 f  = fnunw(f,360.0);
 ms = fnrad(ms);
 mm = fnrad(mm);
 f  = fnrad(f);

 ddjd = (1.734E-1-3.93E-4*t)*sin(ms)
        +2.1E-3*sin(2.0*ms)-4.068E-1*sin(mm)
        +1.61E-2*sin(2.0*mm)-4E-4*sin(3.0*mm)
        +1.04E-2*sin(2.0*f)-5.1E-3*sin(ms+mm)
        -7.4E-3*sin(ms-mm)+4E-4*sin(2.0*f+ms)
        -4E-4*sin(2.0*f-ms)-6E-4*sin(2.0*f+mm)
        +1E-3*sin(2.0*f-mm)+5E-4*sin(ms+2.0*mm);
 
 a1 = fnlif(a);
 b += ddjd+a-a1;
 b1 = fnlif(b);
 a  = a1+b1;
 b  = b-b1;
 
 return getJulianDate()-2415020.0 - a - b;
}

double timestamp::nearestFullMoon()
{
 double  djd0, djd, k, tn, tf, t, t2,b,c,ms,mm,f,ddjd,a1,b1,a,r;
 
 djd0 = julian(Year+1900,1,0,0,0,0);
 djd  = julian(Year+1900,Month+1,Day,0,0,0);
 k    = fnlif((Year+(djd-djd0)/365.0)*12.3685+0.5);
 tn   = k/1236.85;
 tf   = (k+0.5)/1236.85;
 t    = tf;
 k   += 0.5;

 t2 = t*t;a=29.53*k;
 c  = fnrad(166.56+(132.87-9.173E-3*t)*t);
 b  = 5.8868E-4*k+(1.178E-4-1.55E-7*t)*t2
      +3.3E-4*sin(c)+7.5933E-1;
 ms = 359.2242+360.0*fnunw(k/1.236886E1,1.0)
      -(3.33E-5+3.47E-6*t)*t2;
 mm = 306.0253+360.0*fnunw(k/9.330851E-1,1.0)
      +(1.07306E-2+1.236E-5*t)*t2;
 f  = 21.2964+360.0*fnunw(k/9.214926E-1,1.0)
      -(1.6528E-3+2.39E-6*t)*t2;
 ms = fnunw(ms,360.0);
 mm = fnunw(mm,360.0);
 f  = fnunw(f,360.0);
 ms = fnrad(ms);
 mm = fnrad(mm);
 f  = fnrad(f);

 ddjd = (1.734E-1-3.93E-4*t)*sin(ms)
        +2.1E-3*sin(2.0*ms)-4.068E-1*sin(mm)
        +1.61E-2*sin(2.0*mm)-4E-4*sin(3.0*mm)
        +1.04E-2*sin(2.0*f)-5.1E-3*sin(ms+mm)
        -7.4E-3*sin(ms-mm)+4E-4*sin(2.0*f+ms)
        -4E-4*sin(2.0*f-ms)-6E-4*sin(2.0*f+mm)
        +1E-3*sin(2.0*f-mm)+5E-4*sin(ms+2.0*mm);
 
 a1 = fnlif(a);
 b += ddjd+a-a1;
 b1 = fnlif(b);
 a  = a1+b1;
 b  = b-b1;
 
 r = getJulianDate()-2415020.0 - a - b;

 if (r<-29.53)      r += 29.53; 
 if (r<-29.53/2.0)  r += 29.53;
 return r;
}

string timestamp::getNiceTime(const int land)
{
 string str = "";

 if (Hour<10) str = "0";

 str+= toString(Hour);
 str+= ":";

 if (Minute<10) str += "0";

 str+= toString(Minute);

 return str;
}

void  timestamp::setToLastDayOfMonth(int hour, int dayofweek, int month, int year)
{
 int lastDayOfWeek[7] = {1,1,1,1,1,1,1};
 int day = 1;
 bool     finished = false;
 struct   tm TM;
 time_t   dummy;

 while (!finished)
 {
  TM.tm_sec  = 0;
  TM.tm_min  = 0;
  TM.tm_hour = 1;
  TM.tm_mday = day;
  TM.tm_mon  = month-1;
  TM.tm_year = year-1900;
  
  dummy = mktime(&TM);
  TM    = *(localtime(&dummy));

  if (TM.tm_mon ==(int)(month-1)) lastDayOfWeek[TM.tm_wday] = day;
    else finished = true;
   
  day++;
 }
 
 TM.tm_hour = hour;
 TM.tm_mday = lastDayOfWeek[dayofweek];
 TM.tm_mon  = month-1;
 TM.tm_year = year-1900;

 setMembers(TM);
}

int timestamp::firstDayOfYear()
{
 struct tm TM;
 time_t    dummy;
 
 TM.tm_sec  = 0;
 TM.tm_min  = 0;
 TM.tm_hour = 12;
 TM.tm_mday = 1;
 TM.tm_mon  = 0;
 TM.tm_year = Year;  

 dummy = mktime(&TM);
 TM    = *(localtime(&dummy));

 return TM.tm_wday;
} 

int timestamp::getKW()
{ 
 int firstDay = firstDayOfYear();
 if (!firstDay) firstDay = 7;

 return (DayOfYear+firstDay-1-!DayOfWeek)/7+1;
}
  
int timestamp::dayLightSavingTime()
{
 //--- checks whether daylight saving time is currently in use
  
 timestamp begin;
 timestamp end;
 timestamp selectedTime = *this;
 
 //--- day light saving time begins on last sunday in march 3 am
 begin.setToLastDayOfMonth(3,0,3,Year+1900);

 //--- day light saving time ends on last sunday in october 3 am
 end.setToLastDayOfMonth(3,0,10,Year+1900);

 //-- today: begin of daylight saving time
 if (selectedTime.getMonth()==begin.getMonth() && 
     selectedTime.getDay()  ==begin.getDay()) return 2;

 //-- today: end of daylight saving time
 if (selectedTime.getMonth()==end.getMonth() && 
     selectedTime.getDay()  ==end.getDay()) return 3;

 //-- else: just test if daylight saving time is used
 return (selectedTime>=begin && selectedTime<=end);
}

void timestamp::setToEaster()
{
 //--- sets timestamp to easter sunday of year
 //--- e.g. first sunday after 21st of March (1st day of spring)
 //--- AFTER first full moon

 //-- set to 22nd of March (fullmoon on 21st of March doesn´t count!)

 setTimeStamp(0,0,0,22,2,Year);

 bool passedFullMoon = false;

 //-- step forward daywise as long as no full moon has passed
 //-- the time a sunday is meet

 while (!(passedFullMoon && !getDayOfWeek()))
 {
  if (ceil(nearestFullMoon())==0.0) passedFullMoon = true;
  setSeconds((long int) floor(getSeconds()+3600.0*24.0));
 }
}  
