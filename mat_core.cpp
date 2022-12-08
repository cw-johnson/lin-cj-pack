#include "CJVector.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------MATRIX core functions------------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//---------------------------------- Access Operators -----------------------------------//
//Get row
double &CJMatrix::operator[](const int &row)
{
  return data[row * n];
}

const double &CJMatrix::operator[](const int &row) const
{
  return data[row * n];
}
//Get element
double &CJMatrix::operator()(const int &row, const int &col)
{
  return data[calcIndex(row, col)];
}

//Get element, const
const double &CJMatrix::operator()(const int &row, const int &col) const
{
  return data[calcIndex(row, col)];
}
//Get element
double &CJMatrix::getelem(const int &row, const int &col)
{
  return data[calcIndex(row, col)];
}

//Get element, const
const double &CJMatrix::getelem(const int &row, const int &col) const
{
  return data[calcIndex(row, col)];
}
//Caculate element index
int CJMatrix::calcIndex(const int row, const int col) const
{
  if ((row < m) && (col < n))
  {
    if ((mtype == DENSE) || (mtype==SYMMETRIC))
      return row * n + col;
    // else if (mtype == SYMMETRIC)
    // {
    //   int tcol, trow;
    //   if (col > row){
    //     tcol = row;
    //     trow = col;
    //   }
    //   else{
    //     tcol = col;
    //     trow = row;
    //   }
    //   int offset = (trow*(trow+1)) / 2;
    //   return (offset + tcol);
    // }
    else
      return -10;
  }
  else
  {
    cout << "calcIndex Error Row:" << row << " Column: " << col << endl;
    return -11;
  }
}


//Utility Functions
//Print
void CJMatrix::print(){
  if (m>100){
  for (int i = 0; i < 10; i++)
  {
    if (n>3){
      for (int j = 0; j < 3; j++)
        cout << fixed << setprecision(4) << getelem(i,j)<< "      ";
      cout<<"   ...   ";
      for (int j = n-3; j < n; j++)
        cout << fixed << setprecision(4) << getelem(i,j)<< "      ";
      cout << endl;
    }else{
      for (int j = 0; j < n; j++)
        cout << fixed << setprecision(4) << getelem(i,j)<< "      ";
      cout << endl;    
    }
  }
  cout << "              ...          "<<endl; 
  for (int i = m-10; i < m; i++)
  {
    if (n>3){
      for (int j = 0; j < 3; j++)
        cout << fixed << setprecision(4) << getelem(i,j)<< "      ";
      cout<<"   ...   ";
      for (int j = n-3; j < n; j++)
        cout << fixed << setprecision(4) << getelem(i,j)<< "      ";
      cout << endl;
    }else{
      for (int j = 0; j < n; j++)
        cout << fixed << setprecision(4) << getelem(i,j)<< "      ";
      cout << endl;    
    }
  }    
  cout<<endl;
  }
  else{
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      cout << fixed << setprecision(4) << getelem(i,j);
      cout << "      ";
    }
    cout << endl;
  }
  cout << endl;
  }
}


bool dimMatch(const CJMatrix& lhs, const CJMatrix& rhs){
    if((lhs.len == rhs.len) && (lhs.m == rhs.m) && (lhs.n == rhs.n))
        return true;
    else
        return false;
}


//Set All to Value
void CJMatrix::setAll(double val){
    for(size_t i=0; i<len; i++){
        data[i] = val;
    }
}

//Set All to 0
void CJMatrix::zeros(){
    setAll(0.0);
}

//Add random uniform data to all elements of matrix between min and max
void CJMatrix::addRandomData(double min, double max)
{
  uniform_real_distribution<double> myrand(min, max);
  default_random_engine re; // https://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
  for (size_t i = 0; i < len; i++)
  {
    data[i] += myrand(re);
  }
}

