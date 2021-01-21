#ifndef FIND_PEAKS_H  
#define FIND_PEAKS_H  
#include <mpi.h>
#include <iostream>     // cin, cout, endl, cerr
#include <vector>       // vector
#include <map>
#include <limits>       // std::numeric_limits
#include <math.h> 
#include <functional> 
#include <fstream>      // std::ofstream
#include <numeric>      // std::iota
#include <algorithm>
#include <bits/stdc++.h> 

//BOOST
//#ifdef USE_BOOST 
#include <tuple> // for std::tuple and std::make_tuple.
#include <boost/range/adaptors.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/math/special_functions/next.hpp> // For float_distance.
#include <boost/math/special_functions/cbrt.hpp> // For boost::math::cbrt.
#include <boost/math/interpolators/barycentric_rational.hpp>
//#endif 

template <class T>
class FindPeaks 
{
  public : 
    FindPeaks(){};
   ~FindPeaks(){};

    std::vector<int> SetData(T* array, int nArray)
    { 
      // MAXs 
      max_ids = this->FindPeaksIds(array, 0, nArray); 
      //std::cout<<"nMaxs:"<< max_ids.size() <<" \n";

      // MINs 
      std::transform(array, array+nArray, array, [](T &c){return -1*c;});
      min_ids = this->FindPeaksIds(array, 0, nArray);       
      std::transform(array, array+nArray, array, [](T &c){return -1*c;});
      //std::cout<<"nMins:"<< min_ids.size() <<" \n";

      // ALLs  
      max_ids.insert(max_ids.end(), min_ids.begin(), min_ids.end());
      max_ids.push_back(nArray-1); 
      //this->Print(max_ids);

      // increasing order 
      std::sort(max_ids.begin(), max_ids.end(), std::less<int>()); 

      // decreasing order 
      //std::sort(max_ids.begin(), max_ids.end(), std::greater<int>());
      //this->Print(max_ids); 

      //this->GetRanges(array); 

      return max_ids; 
    }; 

    void Print(std::vector<int> v)
    {
      for(int i=0; i<v.size(); i++){ std::cout<< v[i] <<" ";}; std::cout<<"\n"; 
    }

    void Print(std::vector< std::vector<T> > vector) 
    { 
      for(int i=0; i<vector.size(); i++)
      { 
        for(int j=0; j<vector[i].size(); j++) std::cout<< vector[i][j] <<" ";
        std::cout<<"\n";
      }
    }

    std::vector< std::vector<T> > GetRanges(T* array)
    { 
      std::vector< std::vector<T> > Ranges;
      
      for(int i=1; i<max_ids.size(); i++)
      { 
        int b = max_ids[i-1]; // bottom  
        int t = max_ids[i  ]; // top 
        Ranges.push_back( std::vector<T>{array[b], array[t]} );
      }  
      
      return Ranges;
    }


    std::vector< std::vector<int> > GetRanges()
    { 
      std::vector< std::vector<int> > Ranges;
      
      for(int i=1; i<max_ids.size(); i++)
      { 
        int b = max_ids[i-1]; // bottom  
        int t = max_ids[i  ]; // top 
        Ranges.push_back( std::vector<int>{b,t} );
      }  
      
      return Ranges;
    }

  protected : 
    std::vector<int> FindPeaksIds(T *Array1, int start, int end)
    {   
      std::vector<int> IDs; 

      T prevVal = INT_MIN;
      enum{Ascending, Descending} direction = Ascending;
    
      for (int i=start; i<end; i++)
      {   
        T curVal = Array1[i];
        if (prevVal < curVal) // (still) ascending?  
          direction = Ascending;
        else if (prevVal > curVal) // (still) descending? 
        { 
          if (direction != Descending) // starts descending? 
          { 
            IDs.push_back(i-1); 
            direction = Descending;
            //cout << "peak at index " << i-1 << ": " << prevVal << endl;
          }
        }
        // prevVal == curVal is simply ignored...
        prevVal = curVal;
      }
   
     return IDs;
    } 

  private : 
    std::vector<int> max_ids; 
    std::vector<int> min_ids;

}; 


template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v)
{
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}


struct TerminationCondition  {
  bool operator() (double min, double max)  {
    return abs(min - max) <= 0.000001;
  }
};


double InterpolationBoost(double *x_ptr, double *y_ptr, int nrows, double guess)
{
  double root = std::numeric_limits<double>::max();   
//#ifdef USE_BOOST 
  boost::math::barycentric_rational<double> b(x_ptr, y_ptr, nrows);

  using boost::math::tools::bisect;
  double from = x_ptr[0];  // The solution must lie in the interval [from, to], additionally f(from) <= 0 && f(to) >= 0
  double to = x_ptr[nrows-1];

  std::pair<double,double> result = bisect([=](double x) {return b(x);}, from, to, TerminationCondition());
  root = (result.first + result.second) / 2; 

  double abscissa_3 = root ; 
  int inside = (x_ptr[0]<=abscissa_3)&&(abscissa_3<=x_ptr[nrows-1]);  
//#endif  
  return root; 
}


std::vector<double> GetZeros(std::vector<double> X, std::vector<double> Y, bool debug=false)
{
  assert( X.size()==Y.size() );  

  FindPeaks<double> *findPeaks = new FindPeaks<double>();
  std::vector<int> IDs = findPeaks->SetData(Y.data(), Y.size());

  std::vector<double> Roots;
  std::vector< std::vector<int> > Ranges = findPeaks->GetRanges();

  debug = true; 
  if(debug)
  {
    std::cout<<"\t[GetZeros] Peaks found:"<< Ranges.size() <<" \n"; 
    //for(int i=0; i<Ranges.size(); i++) std::cout<<"\t  "<< i <<") "<< Ranges[i][1] - Ranges[i][0] <<" \n";     
    //std::cout<<" \n"; 
  }

  for(int i=0; i<Ranges.size(); i++)
  {
    int b = Ranges[i][0];
    int t = Ranges[i][1];

    std::vector<double> x(X.data()+b, X.data()+t); 
    std::vector<double> y(Y.data()+b, Y.data()+t); 

    if(x.size() > 1) 
    { 
      std::cout<<"\n ------------------------------------- \n";
      int positives = 0;
      int negatives = 0; 
      for(int i=0; i<x.size(); i++)  
      {
        if(y[i]>=0) positives += 1;   
        if(y[i]< 0) negatives += 1;      
      }

      if( (positives==0)||(negatives==0) ) 	
      {
        std::cout<<" no roots...  \n ";
      }
      else
      { 
        std::cout<<"    roots...  \n ";
        if( (positives>1)&&(negatives>1) )  
        {
          double average = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
          //std::cout<<"\nAverage:"<< average <<" \n"; 
          double root = InterpolationBoost(x.data(), y.data(), x.size(), average);
          //std::cout<<"\nRoot:"<< root <<" \n"; 
          Roots.push_back(root);
        }
	else
        {
/* 
for(int i=0, p=0, n=0; i<x.size(); i++)
{
  if(y[i]>=0) p += 1;
  if(y[i]< 0) n += 1;
 
  std::cout<< i <<" "<< x[i] <<" "<< y[i] <<" "; 
  std::cout<< p <<" "; 
  std::cout<< n <<" ";
  std::cout<<"\n";
}
*/
	}
      } // (positives==0)||(negatives==0)  
      std::cout<<"\n ------------------------------------- \n";
    } // if x.size() > 1  
  } // for

  return Roots;
}




/*
EXAMPLE : 
int main()
{
  std::ofstream myfile;
  std::string fname; 

  int N  = 1000; 

  std::vector<double> X(N,std::numeric_limits<double>::max()); 
  for(int i=0;i<N;i++) X[i] = 2*M_PI/(N-1)*i;   

  std::vector<double> Y(N,std::numeric_limits<double>::max());
  for(int i=0;i<N;i++) Y[i] = sin( 2*X[i] ) + cos( 11*X[i] ) ; 

  //for(int i=0;i<N;i++) std::cout<< X[i] <<" "<< Y[i] <<" \n";
  fname="findPeaks.dat";            
  myfile.open( fname.c_str() );
  for(int i=0;i<X.size();i++) myfile<< X[i] <<" "<< Y[i] <<" \n";
  myfile.close();

  FindPeaks<double> *findPeaks = new FindPeaks<double>();   
  std::vector<int> IDs = findPeaks->SetData(Y.data(), Y.size()); 

  fname="foundPeaks.dat"; 
  myfile.open( fname.c_str() );
  for(int i=0;i<IDs.size();i++) myfile<< X[IDs[i]] <<" "<< Y[IDs[i]] <<" \n";
  myfile.close();

  std::vector<double> Roots; 
  std::vector< std::vector<int> > Ranges = findPeaks->GetRanges();
  for(int i=0; i<Ranges.size(); i++) 
  {
    int b = Ranges[i][0]; 
    int t = Ranges[i][1]; 

    std::vector<double> x(X.data()+b, X.data()+t);
    std::vector<double> y(Y.data()+b, Y.data()+t);

    int sing = 0;   
    for(int i=0; i<x.size(); i++) sing += (y[i]>=0)?(1):(-1);

    for(int i=0; i<x.size(); i++) std::cout<< i <<" "<< x[i] <<" "<< y[i] <<" "<< signbit(-y[i]) <<" \n";
    std::cout<<"\n\n\n";

    if( x.size()==abs(sing) ) // all have the same sig 
    {
    } 
    else
    {
      double average = std::accumulate(x.begin(), x.end(), 0.0) / x.size(); 
      //std::cout<<"\nAverage:"<< average <<" \n"; 
      double root = InterpolationBoost(x.data(), y.data(), x.size(), average);  
      //std::cout<<"\nRoot:"<< root <<" \n"; 
      Roots.push_back(root); 
    }
  }

  boost::math::barycentric_rational<double> b(X.data(), Y.data(), X.size());
 
  fname="roots.dat";
  myfile.open( fname.c_str() );
  for(int i=0;i<Roots.size();i++) myfile<< Roots[i] <<" "<< b(Roots[i]) <<" \n";
  myfile.close();


}
*/

/*
FROM : 
  ../E04_BOOST/ 

EXECUTION :
  mpicxx -std=c++0x -I/Users/poderozita/z2019_1/REPOSITORY/BOOST170_1 findPeaks.cxx  && ./a.out > a.dat
  gnuplot -e "set grid; plot [0:2*pi*1.01] 'findPeaks.dat' w l, 'foundPeaks.dat' w p, 'roots.dat' w p; pause -1"

SEE : 
  https://www.tutorialspoint.com/cplusplus-program-to-find-the-peak-element-of-an-array-using-binary-search-approach
  https://www.geeksforgeeks.org/find-a-peak-in-a-given-array/
*/

#endif
// J. MIGUEL ZAVALA AKE. 2019AUG16. STOCKHOLM, SWEDEN. 
