#include "loader_alya_wrapper.h"
#include "loader_alya.hpp"
#include <iostream> // cin, cout, endl, cerr
#include <assert.h>
#include <string>
#include <vector>
#include <stdlib.h>

static loader_alya*  loader_ptr = NULL; 

//=======================================================================||===//
//=======================================================================||===//
extern "C" 
{
  void*
  c_loader_alya_create_array()
  {
    //return  new loader_alya; 
    loader_ptr = new loader_alya(); 
    return loader_ptr; 
  }


  void 
  c_loader_alya_delete() //(void* ptr_class)
  {
    //delete static_cast<loader_alya*>(ptr_class);
    assert(loader_ptr != NULL); 
    //delete loader_ptr; 
    loader_ptr = NULL; 
  } 


  void 
  c_loader_alya_set_ple_data(char* filename, int idx)
  {
    //(static_cast<loader_alya*>(ptr_class) + id)->set_ple_data(filename, idx); 
    //if(ptr_class!=NULL) static_cast<loader_alya*>(ptr_class)->set_ple_data(filename, idx); 

    assert(loader_ptr != NULL); 
    std::string str( string_f2c(filename, 256) ); 
    loader_ptr->set_ple_data(str.c_str(), idx); 
  }

}


//=======================================================================||===//
extern "C" 
{
  int stripperdLength(char* strng, int length)
  {
    int i; 
    i = length-1; 
    while( ((strng[i]==' ')||(strng[i]==0)) && (i>=0) ) strng[i--]='?'; 
    
    return i+1; 
  }

  const char* string_f2c(char* strng, int length)
  {
    int i = length-1; 
    while( ((strng[i]==' ')||(strng[i]==0)) && (i>=0) ) i--; 
    i+=1; 
    
    std::string str; 
    for(int j=0; j<i; j++) str.push_back(strng[j]); 

    return str.c_str(); 
  }


  /*
  // fortran 
  real(c_double)  i1, j1, k1 
  common /ijk/ i1, j1, k1 

  // c/c++
  struct
  {
    double ii, jj, kk; 
  } ijk; 


  void 
  cosa()
  {
    std::cout<<"C_ijk: "; 
    std::cout<< ijk.ii <<" "; 
    std::cout<< ijk.jj <<" "; 
    std::cout<< ijk.kk <<" "; 
    std::cout<<"\n"; 
  }
  */
}

//=======================================================================||===//
//=======================================================================||===//
