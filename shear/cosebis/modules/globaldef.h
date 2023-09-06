#ifndef GLOBALDEF_H
#define GLOBALDEF_H

//Based on original code by Patrick Simon

/// floating data type to be used throughout (float or double)
//  float is default, some programs won't work with double
#define number  double
#define cnumber complex<number>
#define SIZE(a) (int)a.size()


#include <math.h>
#include <complex>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <sys/stat.h>
#include <vector>

#include "tolower.h"

using namespace std;


// pointer to simple function
typedef number (*FUNC)(number);
typedef number (&FUNCREF)(number);

// PI
const number pi             = 4.0*atan(1.0);

// angular units in [RAD]
const number degree         = pi/180.0;
const number arcsec         = degree/3600.0;
const number arcmin         = degree/60.0;

// time units in [s]
const number second         = 1.0;
const number minute         = 60.0;
const number hour           = 60.0*minute;
const number day            = 24.0*hour;
const number week           = 7.0*day;
const number year           = 52.0*week;
const number century        = 100.0*year;

// some cosmology stuff
const number Msun           = 1.989E30; // [kg/solar mass]
const number Gnewton        = 6.6726E-11; //[m^3/kg/s^2]
const number parsec         = 3.085678E16; //[m]
const number Mpc            = parsec*1E6; //[m]
const number H0             = 100.0; //[km/s/Mpc] 
const number rhoCrit        = 3./8./pi*1E3*H0/Mpc*1E3*H0/Mpc/Gnewton; //[kg/m^3]
const number c              = 300000.0;//[km/s]
const number Dh             = c/H0;//[Mpc]
const number hubbleVolume   = Dh*Dh*Dh;//[Mpc^3]
const number hubbleDistance = Dh;//[Mpc]
const number hubbleTime     = 1.0/H0*parsec*1E3;//[s]

/// directory for data tables of functions
const char tablePath[] = "";

/// integration intervals (logscale)
const number lowFreq  = -5.;
const number highFreq = +6.;

/// write number to binary data file
void   writeValue(number,FILE*);

///returns true if the file(s) exists
bool CheckFilesExist(const string& filename);
bool CheckFilesExist(const vector<string> filename_vec);
bool CheckFolderExist(const string& folderName);
int find_closest_index(vector<number> input_vec,number findme);

/// read number from binary data file
double read(FILE*);

/// which parameters are in which ASCII data file column? requires special format in first line
void   dataFileHeader(const char*,const char*,int&,int&);

/// greeting message
void   greeting(const char*,const char*);

/** Formats a string

 In case the length of the string is less than the default size in the
 second argument, the remaining space is filled up with dots.
*/
string formatted(string,int);

void readformatted(const char*,string&);
void readformatted(const char*,number&);
void readformatted(const char*,int&);

///prints out message and exits if some condition is fullfilled
void alarm(bool condition,const char* message, string="");

#endif
