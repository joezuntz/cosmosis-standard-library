#ifndef MATRIX
#define MATRIX

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>

#include "globaldef.h"

 
 //---------------- defines "numbers" and basic operations thereof --------------
 //----- unfortunately not implemented consequently in the code   ---------------- 
 
 #define _zero          0.0
 #define _one           1.0
 #define _half          0.5
 #define _add(a,b)	a+b
 #define _sub(a,b)	a-b
 #define _mul(a,b)	a*b
 #define _div(a,b)	a/b
 #define _sqr(a)        a*a
 #define _sqrt(a)       sqrt(a)
 #define _inverse(a)    1.0/a
 #define _2num(a)       1.0*a
 #define _num2int(a)    (int) (a)

 //--- resolution of imitated functions
 
 #define RESOLUTION 256
   
 number _cos(const number);
 number _sin(const number);
 number _tan(const number);
 number _acos(const number);
 number _asin(const number);
 number _atan(const number);

 int alphaToUnsigned(number);
 int betaToUnsigned(number);
 
 int rotate(int alpha, int beta,
                 int euler1, int euler2, int euler3,
                 bool=false);
                 

/** General class for matrix operations, reduced to functions that do not rely on NR
 * and modified from the original code by Patrick Simon. Modifications by Marika Asgari.
*/
class matrix
 {
  public:
   
   matrix();

   /// constructs 2dim matrix
   matrix(int,int);

   /// constructs vector=one column matrix
   matrix(int);

   matrix(const matrix&);
   ~matrix();

   /// restores matrix in binary file on construction
   matrix(const char*);

   /// frees occupied memory, makes zero size matrix
   void clear();

   /// resize matrix, all present entries get lost
   void resize(int,int);

   /// resize matrix, one column (vector), all present entries get lost
   void resize(int);
      
   /// load matrix with numbers, input array has to be vector of size rows*columns!
   void  load(number*);

   /// loads particular matrix element with new value
   void  load(int,int,number);

   /// loads pparticular vector(=1 column matrix) with element
   void  load(int,number);
    
   /// get value of element in matrix
   number get(int,int) const;

   /// get value of element in vector
   number get(int) const;

   /// exchange matrix by its transpose
   void   transpose();
   
   /// diagonalise matrix (exchange), return transformation matrix: columns=eigenvectors
   ///this is wrong
   ///only works for symmetric matrices
   matrix diagonal();

   /// returns symmetric part of matrix=1/2*(A+A^T)
   matrix symmetric();
   
   /// returns asymmetric part of matrix=1/2*(A-A^T)
   matrix asymmetric();


   /// returns transposed matrix
   matrix t();
   
   /// returns modulus of matrix
   number mod();  


   /// returns maximum element of matrix
   number max();

   /// returns minimum element of matrix
   number min();

   /// returns the condition (=max(w_i)/min(w_i), w_i singular values) of a matrix
   /// is infinte if matrix is singular
   number cond();

   /// returns the trace of the matrix
   number trace();

   /// remove everthing beyond a strip about the diagonal
   void   strip(int);

   /// makes neutral matrix (diagonal=1, remainder zero)
   void    neutral();

   /// loads all elements with same value, by default zero
   void    zero(number=_zero);

   /// formatted print out of matrix to stream
   void    printOut(ostream&,int=4);
   
    /// formatted print out of matrix to stream without a new line at the end
   void    printOutNoNewLine(ostream&,int=4);

   /// formatted print out of matrix to ASCII file
   void    printOut(const char*,int=4);

   /// formatted print out of matrix to ASCII file without a new line at the end
   void    printOutNoNewLine(const char*,int=4);

   /// store matrix in binary format
   void    store(const char*);
   

   /// restore matrix content from binary file, present entries get lost
   void    restore(const char*);
   void    restore(string str) {restore(str.c_str());}

   ///returns sum of all ellements
   number sum();


   /// add a column to matrix from ASCII file
   void    addColumnFromFile(const char*,int,int);

   /// get pointer to element array of matrix
   number* getElements();
   
   /// matrix error message plus abort
   void error(const char*);

   /// returns either element of 1x1 matrix or modulus for larger matrices
   number scalar();   

   /// size of matrix=rows*columns
   int size() const;

   /// get particular row of matrix as vector
   matrix getRow(int);

   /// get particular column of matrix as vector
   matrix getColumn(int);
   

   /// flips matrix about x-axis
   void  flipx();

   /// flips matrix about y-axis
   void  flipy();
   

   /// reads matrix from ASCII-file
   matrix& readFromASCII(const char* filename,int cols=1000, int rows=1000);

   ///this one doens't need to know the cols and rows
   matrix& readFromASCII_marika(const char* filename);
	
	///marika's sub Matix Keep routine, takes the elements which rows and columns should be kept and return the reduced matrix, the input matrix should be square
	matrix subMatrixKeep(vector<int> elements);
	
	///marika's routine does the scalar product of two matrices with the same size
	matrix ScalarProduct(const matrix& mat2);

	///marika's routine, takes the diagonal elements of the square matrix and saves them in a vector
	matrix diag();
	
	///marika's routine, returns matrix with elements equal to the elements of the object matrix to the power of pow
	matrix power(number pow);

  // void    div(const matrix&);    
   bool    equals(const matrix&);
   void    operator=(number*);
      
   void    add(const matrix&);
   void    add(number);
   void    add(int,int,number);
   void    add(int,number);
   void    sub(const matrix&);
   void	  sub(number);
   void    sub(int,int,number);
   void    sub(int,number);
   void    mul(const matrix&);
   void    mul(number);
   void    mul(int,int,number);
   void    mul(int,number);
   void    div(number);
   void    div(int,int,number);
   void    div(int,number);
   // void    sqrt();
   
   void  operator=(const matrix&);      
   void  operator+=(const matrix&);
   void  operator-=(const matrix&);
   void  operator+=(number);
   void  operator*=(const matrix&);
   void  operator*=(number);
   void  operator/=(const matrix&);
   void  operator/=(number);
         
   int columns;
   int rows;   
   int inArray(int,int) const;
   number*  elements;

 };
   
/// add two matrices
matrix operator+(const matrix&, const matrix&);
/// add a constant to every matrix element
matrix operator+(const matrix&,number);
/// add a constant to every matrix element
matrix operator+(number,const matrix&);
 
/// subtract one matrix from another
matrix operator-(const matrix&, const matrix&);
/// subtract a constant from every matrix element
matrix operator-(const matrix&, number);
/// subtract a constant from every matrix element
matrix operator-(number,const matrix&);
 
/// product of two matrices
matrix operator*(const matrix&, const matrix&);
/// multiply every matrix element by a number
matrix operator*(const matrix&, number);
/// multiply every matrix element by a number
matrix operator*(number, const matrix&);
 

/// divide every matrix element by a number
matrix operator/(const matrix&,number);

/// multiply a matrix N times by itself
matrix operator^(const matrix&,int);

/// 3-element vector product!
matrix operator^(const matrix&,const matrix&);

/// compares two matrices 
bool   operator==(const matrix&,const matrix&);

/// feed matrix content formatted to a stream
ostream& operator<<(ostream&, matrix&);

/// 3x3 rotation matrix: about x-axis
matrix D1(number);
/// 3x3 rotation matrix: about y-axis
matrix D2(number);
/// 3x3 rotation matrix: about z-axis
matrix D3(number);
/// 3x3 rotation matrix defined by Euler angles
matrix euler(number,number,number);
/// 3x3 rotation matrix plus scaling defined by vector (vector product)
matrix rot3d(const matrix&);

/// computes the square root of a square matrix, if its real eigenvalues are all positive (not checked!)
//matrix sqrt(const matrix&);

/// converts cartesian to polar coordinates
void polarCoords(number,number,number,number&,number&,number&);

///marika's sub Matix routine, takes the elements which rows and columns should be substracted and return the reduced matrix, the input matrix should be square
matrix subMatrix(matrix input,vector<int> elements);

///marika's sub Matix routine, takes the rows_to_remove and columns_to_remove and removes the rows and columns
//that should be substracted and return the reduced matrix, the input matrix can be any shape
matrix subMatrix_removeRows_and_Columns(matrix input,vector<int> rows_to_remove, vector<int> columns_to_remove);

///marika's routine finds out if i is not equal to any of the elements
bool notEqual(int i,vector<int> elements);

///marika's routine finds out if i is equal to any of the elements
bool Equal(int i,vector<int> elements);


#endif
