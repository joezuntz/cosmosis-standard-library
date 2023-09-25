#include "globaldef.h"

///checks if a file exists. exits if it doesn't
bool CheckFilesExist(const string& filename)
{
  struct stat buffer;   
    return (stat (filename.c_str(), &buffer) == 0); 
}

bool CheckFilesExist(const vector<string> filename_vec)
{
  struct stat buffer;   
  for(int i=0; i<filename_vec.size(); i++)
  {
    if (stat (filename_vec[i].c_str(), &buffer)!=0)
    {
      return false;
    }
  }
  return true;
}

bool CheckFolderExist(const string& folderName)
{
  struct stat sb;
  if( stat(folderName.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
  {
    clog<<"folder: "<<folderName<<" exists."<<endl;
    return true;
  }
  else
  {
      clog<<"folder: "<<folderName<<" does not exist."<<endl;
    return false;
  }
}

//assumes that the vector is sorted
int find_closest_index(vector<number> input_vec,number findme)
{
  if (findme<input_vec[0])
    return 0;
  if (findme>input_vec[input_vec.size()-1])
    return input_vec.size()-1;

  int index;
  for(index=0; ((index<input_vec.size()) && (input_vec[index]<findme)); index++);///finds the position of the galaxy in the sorted vector with z>=z_min
  if((findme-input_vec[index-1])<(input_vec[index]-findme))
    return index-1;
  else
    return index;
}

void writeValue(number value,FILE* handler)
{
  float dummy = value;
  fwrite((void*) &dummy,sizeof(float),1,handler);
}

double read(FILE* fhandler)
{
  float dummy;
  if (!feof(fhandler))
    fread((void*) &dummy,sizeof(float),1,fhandler);
    
  return dummy;
}

void dataFileHeader(const char* filename,const char* columnName,int& column,int& columns)
{
  string cname = string(columnName);
  toLower(cname);

  ifstream fhandler;
  fhandler.open(filename);
  alarm(fhandler.fail(),"datafilehandler: could not open file",filename);

  string format;
  fhandler >> format;
  fhandler.close();

  if (format.find(":",0)==format.npos)
    {
      column = columns = 1;
      clog << "datafilehandler: WARNING! input file "
	   << filename
	   << " does not contain file header?"
	   << endl;

      return;
    }

  columns=0;
  int    pos=0;
  int    oldPos;
  string name;
  bool   found=false;
  do
    {
      columns++;
      oldPos        = pos;
      pos           = format.find(":",pos);
      name          = format.substr(oldPos,pos-oldPos);
      toLower(name);
      if (name==cname)
	{
	  column = columns;
	  found  = true;
	}
      if (pos!=(int)format.npos) pos++;
    } while (pos!=(int)format.npos);

  alarm(!found, 
	(string("datafilehandler: could not find column ")+columnName+
	string(" in file ")+filename).c_str());
}

string formatted(string str, int width)
{
  int n;

  string result=str;

  for(n=str.size();n<=width;n++)
    result+=".";

  result+=" ";

  return result;
}

void  greeting(const char* name,const char* version)
{
    int n;

    clog << endl << endl;

    for(n=0;n<40;n++) clog << "-";
    clog << endl;
    clog << name << endl << endl;
    clog << "contact: " << endl;
    clog << "Patrick Simon <psimon@roe.ac.uk>" << endl;
    clog << "(or psimon@astro.uni-bonn.de)"    << endl << endl;
    clog << "Version " << version << " (compiled)" << endl;
    for(n=0;n<40;n++) clog << "-";
    
    clog << endl << endl;
}

void readformatted(const char* message, string& var)
{
  clog << formatted(message,60);
  cin  >> var;
}

void readformatted(const char* message, number& var)
{
  clog << formatted(message,60);
  cin  >> var;
}

void readformatted(const char* message, int& var)
{
  clog << formatted(message,60);
  cin  >> var;
}

void alarm(bool condition,const char* message, string extra)
{
  if (condition)
    {
      cerr << endl << ">>>>>>>>>> ERROR: " << message << " " << extra << endl << endl;

      exit(1);
    }
}
