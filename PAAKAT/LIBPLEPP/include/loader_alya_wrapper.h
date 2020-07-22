//=======================================================================||===//
//=======================================================================||===//
#ifndef LOADER_ALYA_WRAPPER_H
#define LOADER_ALYA_WRAPPER_H

  typedef struct loader_alya LOADER_ALYA; 

#ifdef __cplusplus
extern "C" // fortran visibility!! 
{ 
#endif

  // loaders arrays
  void*
  c_loader_alya_create_array(); 

  void 
  c_loader_alya_delete(); 

  void 
  c_loader_alya_set_ple_data(char*, int); 

  int stripperdLength(char*, int); 
  const char* string_f2c(char*, int);

  /*
  struct ijk;
  void cosa(); 
  */
  
#ifdef __cplusplus
}
#endif


#endif //LOADER_ALYA_WRAPPER_H
//=======================================================================||===//
//=======================================================================||===//
