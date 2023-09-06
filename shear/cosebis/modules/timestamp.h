#ifndef TIMESTAMP_H
#define TIMESTAMP_H

using namespace std;

#include <string>

double julian(double year,double month, double day, double hour, double minute, double second);
double sgn  (double);
double enlif(double);
double fnunw(double);
double fnrad(double);
double fndeg(double);

/** General class for handling time information

 Extension of the original version of Frank P.
 By Patrick Simon
*/
class timestamp{
    public:
      timestamp(); // sets timestamp to current system-time
      timestamp(const timestamp&);
      timestamp(const int); // Expects seconds defined following the time_t-convention
      timestamp(int, int, int, int, int, int);  // Broken down timestamp
      timestamp(const string& timestampstr); // Expects a string-representation
                                             // as "YYYYMMDDhhmmss"
      void setTimeStamp(int second, int minute, int hour, int day, int month, int year);
      void setTimeStamp(const timestamp&);
      void setSeconds(const long int);
      void parseString(const string&);
      int getSeconds() const;
      string getAsString() const;
      string getFormatted() const;
      void setCurrent();      
      int calculateTMSeconds();
      string getNiceForm(const int land=0) const; 
      string getNiceTime(const int land=0);
      string getMonthStr(int land=0);
      void setTimeStampString();
      int getDay(){ return Day; }
      int getMonth() { return Month; }
      int getYear() { return Year; }
      int getHour() { return Hour; }
      int getMinutes() { return Minute; }
      int getSeconds() { return Seconds; }
      int getSecond()  { return Second; }
      int getDayOfWeek() {return DayOfWeek;}
      int getDayOfYear() {return DayOfYear;}

      //--- Patrick´s implementations
      
      double   getJulianDate();
      double   nearestNewMoon();      	            
      double   nearestFullMoon();
      void     setToLastDayOfMonth(int hour, int dayofweek, int month, int year);
      int      getKW();
      int      firstDayOfYear();
      int      dayLightSavingTime();
      void     setToEaster();
      
      static timestamp getCurrentTime();
      static int toSeconds(const string&);
     
      timestamp& operator=(const timestamp& stamp);

      static const int ENGLISH=0;

    private:
      int Second;
      int Minute;
      int Hour;
      int Day;
      int Month;
      int Year;
      int Seconds;
      int DayOfWeek;
      int DayOfYear;
      
      string TimestampString;      
      
      void setMembers(const struct tm&);
  };

  bool operator!=(const timestamp& t1, const timestamp& t2);
  bool operator==(const timestamp& t1, const timestamp& t2);
  bool operator<(const timestamp& t1, const timestamp& t2);
  bool operator<=(const timestamp& t1, const timestamp& t2);
  bool operator>(const timestamp& t1, const timestamp& t2);
  bool operator>=(const timestamp& t1, const timestamp& t2);
  timestamp operator+(const timestamp& rhs, const int seconds);
  
  
#endif
