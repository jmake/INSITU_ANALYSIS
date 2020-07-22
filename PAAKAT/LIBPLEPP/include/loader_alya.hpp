//#include <vector>
#ifndef LOADER_ALYA_H
#define LOADER_ALYA_H

class loader_alya
{
  public: 
    loader_alya(); 
   ~loader_alya();
    
    void init();
    void end(); 
    
    void get_data(const char*, int=-1); 
    void print_data(); 
    void draw(); 
    void draw_quads(); 
    
    void set_ple_data(const char*, int=-1); 


    // ple 
    int n_vertices; 
    int n_elements; 

    int*    vertex_num_ptr; 
    double* vertex_coords_ptr;

    void get_n_vertices(); 
    void get_n_elements(); 
    void get_vertex_num_ptr();
    void get_vertex_coords_ptr(); 


    //opengl   
    //vector< vector<createpoint> > gl_model; 
    //vector<int>    vertex_num; 
    //vector<double> vertex_coords;
    //std::vector< std::vector<double> > pts; 
}; 

//  #include "loader_alya.cpp"

#endif
