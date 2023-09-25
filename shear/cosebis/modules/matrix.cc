#include "unistd.h"
#include "matrix.h"

 //--------- trigonometric functions

 number _cos(const number a) {return cos(a);}

 number _sin(const number a) {return sin(a);}

 number _tan(const number a) {return tan(a);}

 number _acos(const number a){return acos(a);}

 number _asin(const number a){return asin(a);}

 number _atan(const number a){return atan(a);}

 //----------------------------------------------------------------------
void matrix_error(const char error_text[])
{
  cerr << "Error"<<endl;
  cerr << error_text << endl;
  cerr << ", Closing now\n"<<endl;
  exit(1);
}

 int alphaToUnsigned(number x)
 {
  return (int) (x*RESOLUTION/(2.0*pi));
 }

int betaToUnsigned(number x)
{
  return (int) ((x+pi/2.0)*RESOLUTION/pi);
}
 
int rotate(int alpha,  int beta,
                 int euler1, int euler2, int euler3, 
                 bool deinstall)
{
  static int *lookuptable=NULL;
  static bool wasHere = false;

  if (deinstall)
  {
   //--- use deinstall=true to free allocated memory !

   if (lookuptable==NULL) delete[] lookuptable;

   wasHere = false;

   cerr << "deinitialize lookuptable: rotate" << endl;

   return 0;
  }

  if (wasHere)
  {
   int alpha1, beta1, ab2;
   int alpha2, beta2, ab3;
   int alpha3, beta3;

   //--- use buffer to calculate rotation fast

   //-- convert to int according to resolution
  
   alpha1 = (alpha+euler1)%RESOLUTION;
   beta1  = beta;

   ab2    = lookuptable[alpha1+beta1*RESOLUTION];
   alpha2 = ab2&(RESOLUTION-1);
   beta2  = ab2-alpha2;
   alpha2 = (alpha2+euler2)%RESOLUTION;

   ab3    = lookuptable[alpha2+beta2*RESOLUTION];
   alpha3 = ab3&(RESOLUTION-1);
   beta3  = ab3-alpha3;
   alpha3 = (alpha3+euler3)%RESOLUTION;

   //--- returns direction (alpha,beta) as compact int
   //--- = alpha+beta*RESOLUTION
   //--- take care to convert your texture to this format and you
   //--- can use this as texture index right away !

   return beta3+alpha3;
  }

  clog << "initializing lookuptable: rotate" << endl;
  
  //--- initialize buffer on first call

  lookuptable = new int[RESOLUTION*RESOLUTION];

  int alpha_int;
  int beta_int;
  int alpha_int2;
  int beta_int2;
  number    alpha1;
  number    alpha2;
  number    beta1;
  number    beta2;
  number    x,y,z;

  for(alpha_int=0;alpha_int<RESOLUTION;alpha_int++)
  for(beta_int =0;beta_int <RESOLUTION;beta_int++)
  {
   //--- direction in spherical rads

   alpha1 = alpha_int*(2.0*pi)/RESOLUTION;
   beta1  = beta_int*pi/RESOLUTION-pi/2.0;

   //--- direction in cartesian coords (rotated about x-axis)

   x      = cos(alpha1)*cos(beta1);
   y      = sin(beta1);
   z      = sin(alpha1)*cos(beta1);

   //--- spherical rads in rotated coordinate system
    
   beta2  = asin(z);
   alpha2 = acos(x/cos(beta2));
   if (y/cos(beta2)<0) alpha2 = 2.0*pi-alpha2;

   //--- and converted to integer according to resolution

   alpha_int2 = alphaToUnsigned(alpha2);
   beta_int2  = betaToUnsigned(beta2);

   //--- write to "transfer"-matrix
 
   lookuptable[alpha_int+beta_int*RESOLUTION] = 
    alpha_int2+beta_int2*RESOLUTION;     
  }   
  
  wasHere = true;
  
  return 0;
}

 /*
 
   ----------------------------------------------------------------------
   ------ class matrix --------------------------------------------------
   ----------------------------------------------------------------------
  
*/

ostream& operator<<(ostream& s, matrix& a)
{
  a.printOut(s);
 
  return s;
} 

matrix operator+(const matrix& a, number n)
{
  matrix c=a;
 
  c.add(n);

  return c;
}

matrix operator+(number n, const matrix& a)
{
  matrix c=a;
 
  c.add(n);

  return c;
 }

matrix operator-(const matrix& a, number n)
{
  matrix c=a;
 
  c.sub(n);

  return c;
}

matrix operator-(number n, const matrix& a)
{
  matrix c=a;
 
  c.sub(n);

  return c;
}

matrix operator^(const matrix& a,const matrix& b)
{
  matrix c=rot3d(a);

  c.mul(b);

  return c;
 }

matrix operator^(const matrix& a,int n)
{
  matrix   c=a;
  int i;

  for(i=1;i<n;i++) c.mul(a);

  return c;
 }  

 matrix operator+(const matrix& a, const matrix& b)
 { 
  matrix c=a;

  c.add(b);

  return c;
 }

 matrix operator-(const matrix& a, const matrix& b)
 {  
  matrix c=a;

  c.sub(b);

  return c;
 }

 matrix operator*(const matrix& a, const matrix& b)
 {
  matrix c=a;

  //--- standard product for matrices

  if (a.columns==b.rows)
  {
   c.mul(b);

   return c;
  }

  //--- yields standard scalar product of two vectors
  //--- and an generalization for matrices of same number of columns/rows

  if (a.columns==b.columns || a.rows==b.rows)
  {
   c.transpose();
   c.mul(b);

   return c;
  }

  //--- returns simply a on failure

  c.error("incompatible matrices: mul");

  return a;
 }

 matrix operator*(const matrix& a, number n)
 {
  matrix c=a;

  c.mul(n);

  return c;
 }

 matrix operator*(number n, const matrix& a)
 {
  matrix c=a;

  c.mul(n);

  return c;
 }


 // matrix operator/(const matrix& a, const matrix& b)
 // {
 //  matrix c=a;
 //  matrix d=b;

 //  c.mul(d.inverse());

 //  return c;
 // }

 // matrix operator/(number n, const matrix& a)
 // {
 //  matrix c=a;
  
 //  c.inverse();
 //  c.mul(n);

 //  return c;
 // }

 matrix operator/(const matrix& a, number n)
 {
  matrix c=a;

  c.div(n);

  return c;
 }

 bool operator==(matrix& a, const matrix& b)
 {
  return a.equals(b);
 }

 //-------- matrix methods ----------------------------------------------

 number* matrix::getElements() {return elements;}

 number matrix::scalar()
 {
  if (columns==1 && rows==1) return elements[0];
   else return mod();
 }

 void matrix::error(const char *message)
 {
  cerr << message << endl;

  exit(1);
 }

 matrix::matrix(int x, int y)
 {
  columns = x;
  rows    = y;

  elements= new number[x*y];
  neutral(); 
 }

matrix::matrix(const char* file)
{
  columns = 0;
  rows    = 0;
  elements= NULL;
  restore(file);
}

void matrix::resize(int x,int y)
{
  clear();

  columns = x;
  rows    = y;

  elements= new number[x*y];

  neutral();
}

void matrix::resize(int x)
{
  resize(1,x);
  neutral();
}

matrix::matrix(int y)
 {
  columns = 1;
  rows    = y;

  elements = new number[y];
  neutral();
 }

matrix::matrix()
 {
  columns = 0;
  rows    = 0;
  
  elements = NULL;
 }

matrix::matrix(const matrix& a)
 {
  columns  = a.columns;
  rows     = a.rows;
  elements = new number[rows*columns];
  
  memcpy(elements,a.elements,sizeof(number)*columns*rows);
 }

 void matrix::operator=(const matrix& a)
 {
  if (elements!=NULL) delete[] elements;

  columns = a.columns;
  rows    = a.rows;
  elements= new number[columns*rows];

  memcpy(elements,a.elements,sizeof(number)*columns*rows);
 }

matrix::~matrix() 
 {
   clear(); 
 }

void matrix::clear()
 {
  if (elements!=NULL) delete[] elements;

  elements = NULL;
  columns  = 0;
  rows     = 0;
 }

 void matrix::load(number *newElements)
 {
   memcpy(elements,newElements,sizeof(number)*columns*rows);
 }



 void matrix::load(int x,int y,number value)
 {
    // if (x>=columns)
    //   error("too many columns in matrix");
    // if (y>=rows)
    //   error("too many rows in matrix");
    elements[inArray(x,y)] = value;
 }

 void matrix::load(int x,number value)
 {
   // //clog<<"value"<<endl;
   // if (x>=columns)
   //    error("too many columns in matrix");
   load(x,0,value);
 }  

 void matrix::operator=(number *newElements)
 {
    load(newElements);
 }

int matrix::inArray(int x,int y) const
{
  if (!x && x>=columns || !y && y>=rows)
  {
    cerr << "matrix: invalid reference to element: "
	    << " x="  << x
	    << " y=" << y
	    << endl;

    x%=columns;
    y%=rows;
  }

  return x+y*columns;
}

number matrix::get(int x) const
{
    return get(x,0);
}

number matrix::get(int x,int y) const
{
  return elements[inArray(x,y)];
}

void matrix::add(const matrix& b)
{
  int n;

  if (b.columns*b.rows==columns*rows)
    for(n=0;n<columns*rows;n++) 
      elements[n]+= b.elements[n];
  else error("incompatible matrices: add");

}

void matrix::add(number n)
{
  matrix c(columns,rows);

  add(c*n);
}

 void matrix::add(int x,int y,number n)
 {
   elements[x+y*columns]+=n; 
 }

 void matrix::add(int x,number n)
 {
   elements[x]+=n;
 }

 void matrix::sub(const matrix& b)
 {
  int n;

 if (b.columns*b.rows==columns*rows)
   for(n=0;n<columns*rows;n++) 
     elements[n]-= b.elements[n];
 else error("incompatible matrices: sub");  

 }

 void matrix::sub(number n)
 {
  matrix c(columns,rows);

  sub(c*n);
 }

 void matrix::sub(int x,int y,number n)
 {
   elements[x+y*columns]-=n;
 }

 void matrix::sub(int x,number n)
 {
   elements[x]-=n;
 }

 void matrix::mul(const matrix& b)
 {
  if (b.rows!=columns)
   error("incompatible matrices: mul");
  
  number*  newElements = new number[b.columns*rows];
  long     col,row,i;
  number   result;

  for(row=0;row<(long)rows;row++)
    for(col=0;col<(long)b.columns;col++)
    {
     for(result=0.,i=0;i<(long)columns;i++)
       result+= 
         elements[i+columns*row]*b.elements[col+b.columns*i];   

     newElements[col+row*b.columns] = result;
    }

  columns = b.columns;
  if (elements!=NULL) 
    delete[] elements;

  elements = newElements;
 }

 void matrix::mul(number f)
 {
  int n;

  for(n=0;n<columns*rows;n++) 
    elements[n]*= f;
 }

 void matrix::mul(int x,int y,number n)
 {
   elements[x+y*columns]*=n;
 }

 void matrix::mul(int x,number n)
 {
   elements[x]*=n;
 }

 void matrix::div(number f)
 {
  int n;

  for(n=0;n<columns*rows;n++) 
    elements[n]/= f;
 }

// void matrix::sqrt()
//  {
//   int n;

//   for(n=0;n<columns*rows;n++) 
//     elements[n]=sqrt(elements[n]);
//  }

 void matrix::div(int x,int y,number n)
 {
   elements[x+y*columns]/=n;
 }

 void matrix::div(int x,number n)
 {
   elements[x]/=n;
 }

 bool matrix::equals(const matrix& b)
 {
  if (columns!=b.columns || rows!=b.rows) return false;

  int n;
  for(n=0;n<columns*rows;n++) 
    if (elements[n]!=b.elements[n]) return false;

  return true;
 }

 void matrix::transpose()
 {
  number *newElements = new number[rows*columns];

  int n;
  int m;

  for(n=0;n<rows;n++)
    for(m=0;m<columns;m++)  
      newElements[n+rows*m] = 
	elements[m+n*columns];

  memcpy(elements,newElements,sizeof(number)*columns*rows);
  delete[] newElements;

  n       = rows;
  rows    = columns;
  columns = n;
 }

 void matrix::zero(number fill)
 {
  int n;

  for(n=0;n<rows*columns;n++) 
    elements[n] = fill;
 }

number matrix::sum()
{
  number result=0.;
  for(int n=0;n<rows*columns;n++) 
    result+=elements[n];
  return result;
}


 void matrix::neutral()
 {
   long m,n;
   for(n=0;n<(long)rows;n++)
     for(m=0;m<(long)columns;m++)
       elements[m+n*columns] = (m==n ? 1. : 0.);
 }  

 matrix matrix::symmetric()
 {
  matrix c(columns,rows);
  
  if (columns!=rows) 
    return c;

  long m,n;
  for(n=0;n<(long)rows;n++)
    for(m=0;m<(long)columns;m++)
      c.elements[m+n*columns] = 
	elements[m+n*columns]+elements[n+m*columns];
  
  c.mul(_half);

  return c;
 }

 matrix matrix::asymmetric()
 {
  matrix c(columns,rows);
  
  if (columns!=rows) return c;

  long m,n;
  for(n=0;n<(long)rows;n++)
    for(m=0;m<(long)columns;m++)
      c.elements[m+n*columns] = 
	elements[m+n*columns]-elements[n+m*columns];

  c.mul(_half);

  return c;
 }
  
 number matrix::mod()
 {
  number result = _zero;

  int n; 
  for(n=0;n<columns*rows;n++) 
    result = _add(result,_sqr(elements[n]));

  return _sqrt(result);
 }


/*
 matrix matrix::norm()
 {
  matrix c(columns,rows);
  c.load(elements);
  c.normalize();

  return c;
 }*/
/*
 void matrix::normalize()
 {
  number modulus = mod();

  if (modulus) 
    div(modulus);
 }*/

number  matrix::trace()
{
  number trace=0.;

  if (rows!=columns)
    return 0.;

  for(int n=0;n<rows;n++)
    trace += get(n,n);

  return trace;
}

void matrix::strip(int d)
{
  int n,m;
  
  for(n=0;n<rows;n++)
    for(m=0;m<columns;m++)
      if (abs(m-n)>d)
	load(m,n,0.);
}

 matrix matrix::t()
 {
  matrix b(columns,rows);
  b.load(elements);

  b.transpose();

  return b;
 }

 number matrix::max()
 {
   number maximum=0.;
   number value=0.;
   int n;
   int m;

   for(n=0;n<rows;n++)
     for(m=0;m<columns;m++)
       if (!n && !m)
	 maximum = get(m,n);
       else if ((value=get(m,n))>maximum) maximum = value;

   return maximum;
 }


 number matrix::min()
 {
   number minimum=0.;
   number value=0.;
   int n;
   int m;

   for(n=0;n<rows;n++)
     for(m=0;m<columns;m++)
       if (!n && !m)
	 minimum = get(m,n);
       else if ((value=get(m,n))<minimum) minimum = value;

   return minimum;
 }

 void  matrix::operator+=(const matrix& a)
 {
  matrix b = *this+a;

  *this=b;
 }

 void  matrix::operator-=(const matrix& a)
 {
  matrix b = *this-a;

  *this=b;
 }
   
 void  matrix::operator+=(number a)
 {
  matrix b = *this+a;
 
  *this=b;
 }

 void  matrix::operator*=(const matrix& a)
 {
  matrix b = *this*a;

  *this = b;
 }

 void  matrix::operator*=(number a)
 {
  matrix b = *this*a;

  *this = b;
 }
 /*
 void  matrix::operator/=(const matrix& a)
 {
  matrix b = *this/a;

  *this = b;
 }

 void  matrix::operator/=(number a)
 {
  matrix b = *this/a;

  *this = b;
 }
*/
void matrix::printOut(const char* filename,int prec)
{
  ofstream fhandler;
  fhandler.open(filename);
  if (!fhandler.fail())
    {
      printOut(fhandler,prec);
      fhandler.close();
    }
  else
    cerr << "matrix: could not open " << filename << " for writing" << endl;
}

void matrix::printOutNoNewLine(const char* filename,int prec)
{
  ofstream fhandler;
  fhandler.open(filename);
  if (!fhandler.fail())
    {
      printOutNoNewLine(fhandler,prec);
      fhandler.close();
    }
  else
    cerr << "matrix: could not open " << filename << " for writing" << endl;
}


void matrix::printOut(ostream& out, int prec)
{
  int n;
  int m;

  out.setf(ios::scientific | ios::showpos);
  out.precision(prec);
  
  for(n=0;n<rows;out << "\n",n++)
  for(m=0;m<columns;m++)
    out << elements[m+n*columns] << "\t";
  out << endl;
}  
/*
//testing
void matrix::printOut(ostream& out, int prec)
{
  int n;
  int m;

  clog.setf(std::ios::scientific | std::ios::showpos);
  clog.precision(prec);
  clog<<"in printOut"<<endl;
  for(n=0;n<rows;out << "\n",n++)
  for(m=0;m<columns;m++)
    clog << elements[m+n*columns] << "\t";
  clog << endl;
}  */

void matrix::printOutNoNewLine(ostream& out, int prec)
{
  int n;
  int m;

  out.setf(ios::scientific | ios::showpos | ios::floatfield);
  out.precision(prec);
  
  for(n=0;n<rows;out << "\n",n++)
  for(m=0;m<columns;m++)
    out << elements[m+n*columns] << "\t";
  //out << endl;
}  


void matrix::store(const char* filename)
{
  FILE* fhandler = fopen(filename,"wb");

  if (fhandler==NULL)
  {
      cerr << "matrix: could not open filename " 
	         << filename<< endl;
      return;
  }

  fwrite((void*) &columns,sizeof(int),1,fhandler);
  fwrite((void*) &rows,sizeof(int),1,fhandler);
  
  number tmp;
  for(int n=0;n<columns*rows;n++)
  {
      tmp = (number) elements[n];
      fwrite((void*) &tmp,sizeof(number),1,fhandler);
  }

  fclose(fhandler);
}

void matrix::restore(const char* filename)
{
  clear();

  FILE* fhandler = fopen(filename,"rb");

  if (fhandler==NULL)
  { 
     clog<<"matrix: not able to restore: "<<filename<<endl;
	  exit(1);
     return;
  }

  fread((void*) &columns,sizeof(int),1,fhandler);
  fread((void*) &rows,sizeof(int),1,fhandler);

  elements = new number[columns*rows];
  number tmp;
  for(int n=0;n<columns*rows;n++)
  {
      fread((void*) &tmp,sizeof(number),1,fhandler);
      elements[n] = (number) tmp;
  }
      
  fclose(fhandler);
}


int matrix::size() const
{
  return columns*rows;
}

 matrix matrix::getRow(int row)
 {
   matrix result(columns,1);
   int n;

   for(n=0;n<columns;n++)
     result.load(n,row,get(n,row));

   return result;
 }

 matrix matrix::getColumn(int column)
 {
   matrix result(1,rows);
   int n;

   if (column>=columns)
     {
       cerr << "matrix: out of area column grab [STOPPED]!"
	    << endl;

       cin >> n;

       return result;
     }

   for(n=0;n<rows;n++)
     result.load(0,n,get(column,n));

   return result;
 }

void matrix::addColumnFromFile(const char* filename,int col,int cols)
{
  ifstream input;
  input.open(filename);
  if (input.fail())
    {
      cerr << "matrix: could not open file "
	   << filename
	   << endl;

      exit(1);
    }

  // skip header
  string dummy;
  input >> dummy;
 
  // determine number of lines of file
  int    n;
  long   lines = 0;
  number x;
  while (!input.eof())
    {
      for(n=0;n<cols;n++)
	input >> x;

      lines+=(!input.eof());
    }
  input.close();

  // set up result matrix
  int newCols = columns+1;
  int newRows = (rows>lines ?rows : lines);
  matrix result(newCols,newRows);
  result.zero();

  int m;
  if (rows&&columns)
    for(n=0;n<rows;n++)
      for(m=0;m<columns;m++)
	result.load(m,n,get(m,n));

  // now add new column from file
  number v=0.;
  ifstream input2;
  input2.open(filename);
  if (input2.fail())
    {
      cerr << "matrix: could not reopen file "
	   << filename
	   << endl;
      
      exit(1);
    }

  // skip header, again...
  input2 >> dummy;

  lines=0;
  while (!input2.eof())
    {
      for(n=1;n<=cols;n++)
	{
	  input2 >> x;
	  if (n==col)
	    v = x;
	}

      if (!input2.eof()) 
	result.load(newCols-1,lines,v);
      lines++;
    }
  input2.close();

  *this=result;
}


void matrix::flipx()
{
  for(int y=0;y<=(int)floor(rows/2-1);y++)
    for(int x=0;x<columns;x++)
    {
    	number tmp = get(x,y);
    	load(x,y,get(x,rows-y-1));
    	load(x,rows-y-1,tmp);
    }
}

void matrix::flipy()
{
  for(int x=0;x<=(int)floor(columns/2-1);x++)
    for(int y=0;y<rows;y++)
      {
	number tmp = get(x,y);
	load(x,y,get(columns-x-1,y));
	load(columns-x-1,y,tmp);
      }
}

matrix matrix::subMatrixKeep(vector<int> elements)
{
	if(columns!=rows)
	{
		clog<<"columns and rows are not equal"<<endl;
		exit(1);
	}
	int subMColumns=elements.size();
	//clog<<"columns="<<input.columns<<"  elements.size()="<<elements.size()<<"  subMColumns="<<subMColumns<<endl;
	matrix subM(subMColumns,subMColumns);

	int m=0;
// 	for(int i=0;i<elements.size();i++)
// 		clog<<"elements["<<i<<"]="<<elements[i]<<endl;
	for(int i=0; i<columns;i++)
	{
		int n=0;
		if(Equal(i,elements))
		{
			for(int j=0;j<columns;j++)
			{
				if(Equal(j,elements))
				{
					subM.load(m,n,get(i,j));
// 					clog<<"m="<<m<<"  n="<<n<<" i="<<i<<" j="<<j<<endl;
					n++;
				}
			}
			m++;
		}
	}
	return subM;
}


matrix subMatrix(matrix input,vector<int> elements)
{
	if(input.columns!=input.rows)
	{
		clog<<"columns and rows are not equal"<<endl;
		exit(1);
	}
	int subMColumns=input.columns-elements.size();

	matrix subM(subMColumns,subMColumns);

	int m=0;
// 	for(int i=0;i<elements.size();i++)
// 		clog<<"elements["<<i<<"]="<<elements[i]<<endl;
	for(int i=0; i<input.columns;i++)
	{
		int n=0;
		if(notEqual(i,elements))
		{
			for(int j=0;j<input.columns;j++)
			{
				if(notEqual(j,elements))
				{
					subM.load(m,n,input.get(i,j));
					n++;
				}
			}
		m++;
		}
	}
	return subM;
}


matrix subMatrix_removeRows_and_Columns(matrix input,vector<int> rows_to_remove, vector<int> columns_to_remove)
{
  int subMRows=input.rows-rows_to_remove.size();
  int subMColumns=input.columns-columns_to_remove.size();
  clog<<"rows="<<input.rows<<"  rows_to_remove.size()="<<rows_to_remove.size()<<"  subMRows="<<subMRows<<endl;
  matrix subM(subMColumns,subMRows);

  int m=0;
  for(int i=0; i<input.columns;i++)
  {
    int n=0;
    if(notEqual(i,columns_to_remove))
    {
      for(int j=0;j<input.rows;j++)
      {
        if(notEqual(j,rows_to_remove))
        {
          subM.load(m,n,input.get(i,j));
          n++;
        }
      }
    m++;
    }
  }
  return subM;
}

bool notEqual(int i,vector<int> elements)
{
	bool result=true;
	for(int j=0;j<elements.size();j++)
		if(i==elements[j])
			result=false;
	return result;
}

bool Equal(int i,vector<int> elements)
{
	bool result=false;
	for(int j=0;j<elements.size();j++)
		if(i==elements[j])
			result=true;
	return result;
}

matrix& matrix::readFromASCII(const char* filename, int cols, int rows)
{
	ifstream fhandler(filename);
	
	alarm(fhandler.fail(),"could not open file",filename);
	if(fhandler.fail())
		exit(1);

	// read in; number of rows may be smaller than rows
	clear();
	vector<matrix> store;
	for(int y=0;y<rows&&!fhandler.fail()&&!fhandler.eof();y++)
	{
		matrix tmpmat(cols);
		for(int x=0;x<cols;x++)
		{	  
			string tmp;
			fhandler >> tmp;
			tmpmat.load(x,atof(tmp.c_str()));
		}

      if (!fhandler.fail()&&!fhandler.eof())
			store.push_back(tmpmat);
	}
	fhandler.close();

	// transfer to matrix array
	resize(cols,store.size());
	for(int y=0;y<store.size();y++)
		for(int x=0;x<cols;x++)
			load(x,y,store[y].get(x));

  return *this;
}

//reads in ascii files. Doesn't need to know how many lines there are.
matrix& matrix::readFromASCII_marika(const char* filename)
{
  ifstream fhandler(filename);
  
  alarm(fhandler.fail(),"could not open file",filename);
  if(fhandler.fail())
    exit(1);

  // read in; number of rows may be smaller than rows
  clear();

  string line;
  int counter=0;
  int cols_temp=0,cols=0;
  // vector<stringstream> stream_vec;
  vector< vector<number> > temp_vec_vec;
  while(getline(fhandler,line))
  {
    stringstream stream(line);
    // string str=string(stream.peek());
    // if (str=="#")
    // {
    //   //nColumns++;
    //   string line;
    //   getline(fhandler,line);
    //   clog<<line<<endl;
    // }
    // else
    // {
      //stream_vec.push_back(stream);
      cols_temp=0;
      number temp;
      vector<number> temp_vec;
      while(stream>>temp)
      {
        temp_vec.push_back(temp);
        //clog<<temp<<endl;
        cols_temp++;
      }
      if(cols_temp>0)
      {
        cols=cols_temp;
        temp_vec_vec.push_back(temp_vec);
        counter++;
      }
      else
      {
        clog<<"empty line"<<endl;
      }
      //clog<<"columns="<<cols<<endl;
      //clog<<line_vec[counter]<<endl;
       
    // }
  }
  fhandler.close();
  rows=counter;
  //cols=columns;
  clog<<"rows="<<rows<<endl;
  clog<<"columns="<<cols<<endl;

  // transfer to matrix array
  resize(cols,rows);
  for(int y=0;y<rows;y++)
  {
    for(int x=0;x<cols;x++)
    {
      load(x,y,temp_vec_vec[y][x]);
      //clog<<get(x,y)<<"\t";
    }
    //clog<<endl;
  }

  return *this;
}


matrix euler(number a,number b,number c)
 {
  return D1(a)*D2(b)*D1(c);
 }

 matrix D1(number angle)
 {
  matrix a(3,3);

  number d1[] = {_cos(angle),-_sin(angle),_zero,
                 _sin(angle),_cos(angle),_zero,
                 _zero,_zero,_one};
  
  a.load(d1);

  return a;
 }

 matrix D2(number angle)
 {
  matrix a(3,3);

  number d2[] = {_cos(angle),_zero,+_sin(angle),
                 _zero,_one,_zero,
                 -_sin(angle),_zero,_cos(angle)};

  a.load(d2);

  return a;
 }

 matrix D3(number angle)
 {
  matrix a(3,3);

  number d3[] = {_one,_zero,_zero,
                 _zero,_cos(angle),-_sin(angle),
                 _zero,_sin(angle),_cos(angle)};
  
  a.load(d3);

  return a;
 }

 matrix rot3d(const matrix& a)
 {
  number plugIn[] = {0,-a.elements[2],+a.elements[1],
                     +a.elements[2],0,-a.elements[0],
                     -a.elements[1],a.elements[0],0};

  matrix c(3,3);
  c.load(plugIn);
 
  return c;
 }


void polarCoords(number x,number y,number z,
                 number& b, number& l,number& r)
{
 number rho = sqrt(x*x+y*y);

 r   = sqrt(x*x+y*y+z*z);

 if ((rho==0.0) && (z>0))    b = 90.0;
  else 
 if ((rho==0.0) && (z==0.0)) b = 0.0;
  else
 if ((rho==0.0) && (z<0))    b = -90.0;
  else                       b = 180.0/pi*_atan(z/rho);

 if ((x==0.0) && (y==0.0))   l = 0.0;
  else
  {
   number phi = 2.0*180.0/pi*_atan(y/(fabs(x)+rho));

   if ((x>=0.0) && (y>=0.0)) l = phi;
    else
   if ((x>=0.0) && (y<0.0))  l = 360.0+phi;
    else                     l = 180.0-phi;
  }
}


matrix matrix::ScalarProduct(const matrix& mat2)
{
  if(rows!=mat2.rows || columns!=mat2.columns){
    matrix_error("dimentions don't match for scalar product");    
  }

	matrix product_mat(columns,rows);

	for(int r=0; r<rows; r++)
		for(int c=0; c<columns; c++)
			product_mat.load(c,r,get(c,r)*mat2.get(c,r));

	return product_mat;
}

matrix matrix::diag()
{
	if(rows!=columns)
		matrix_error("not a square matrix, can't take the diagonal elements");
	matrix vector_mat(rows);
	for(int i=0; i<rows; i++)
		vector_mat.load(i,get(i,i));
	return vector_mat;
}

matrix matrix::power(number power)
{
	matrix power_mat(columns,rows);
	for(int i=0; i<columns; i++)
		for(int j=0; j<rows; j++)
			power_mat.load(i,j,pow(get(i,j),power));
	return power_mat;
}
