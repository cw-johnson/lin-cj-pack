#include "CJVector.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//--------------------------------------VECTOR core functions------------------------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Utility Functions
//Print
void CJVector::print(){
    for(size_t i=0;i<len;i++)
        cout<<data[i]<<"  ";
    cout<<endl;
}

//Output Stream Operator Overload
ostream& operator<<(ostream& os, const CJVector &res){
    for(size_t i=0;i<res.len;i++)
        os<<res.data[i]<<"  ";
    os<<endl;
    return os;
}
//Set All to Value
void CJVector::setAll(double val){
    for(size_t i=0; i<len; i++){
        data[i] = val;
    }
}

//Set All to 0
void CJVector::zeros(){
    setAll(0.0);
}

//Add random uniform data to all elements of matrix between min and max
void CJVector::addRandomData(double min, double max)
{
  uniform_real_distribution<double> myrand(min, max);
  default_random_engine re; // https://stackoverflow.com/questions/2704521/generate-random-double-numbers-in-c
  for (size_t i = 0; i < len; i++)
  {
    data[i] += myrand(re);
  }
}

bool dimMatch(const CJVector& lhs, const CJVector& rhs){
    if(lhs.len == rhs.len)
        return true;
    else
        return false;
}