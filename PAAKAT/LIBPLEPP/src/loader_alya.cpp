// General 
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector 
#include <cstdlib>
#include <cstring>
#include "read_file.hpp"
#include "loader_alya.hpp"
//#include <GL/glew.h> // libglewX.X-dev
using namespace std;

//=======================================================================||===//
class createpoint
{
  public: 
    double x; 
    double y;
    double z;

    void init(vector<double> p, int shift = 0)
    { 
      x = p[0+shift]; 
      y = p[1+shift]; 
      z = p[2+shift]; 
    }; 
}; 
//=======================================================================||===//


//=======================================================================||===//
//=======================================================================||===//
loader_alya::loader_alya() 
{
}


loader_alya::~loader_alya() 
{
}


void
loader_alya::init() 
{
  vertex_num_ptr    = NULL; 
  vertex_coords_ptr = NULL; 
}


void 
loader_alya::end()
{
  if(vertex_num_ptr    != NULL) delete vertex_num_ptr; 
  if(vertex_coords_ptr != NULL) delete vertex_coords_ptr; 
}

void loader_alya::get_data(const char* cname, int Id)
{
    string fdata = string(cname); 
    
    // Reading 
    read_log_file Data[3]; 
    Data[0].set_name(fdata+"_coords.alya"); 
    Data[0].run();
    vector< vector<double> > vcoords( Data[0].get_vdata() );
    Data[0].end(); 

    Data[1].set_name(fdata+"_faces_nodes.alya"); 
    Data[1].run();
    vector< vector<double> > vfaces_nodes( Data[1].get_vdata() );
    Data[1].end(); 

    Data[2].set_name(fdata+"_faces_ids.alya"); 
    Data[2].run();
    vector< vector<double> > vfaces_ids( Data[2].get_vdata() );
    Data[2].end(); 

    // Choosing local mesh 
    int col00 = 1; 
    n_elements = 0;
    vector<int>    vertex_num; 
    vector<double> vertex_coords;

    vector< vector<createpoint> > gl_model; 
    for(int i=0; i<vfaces_ids.size(); i++)
    {
      if(vfaces_ids[i][col00] == Id)
      {
        n_elements++;
        int n_cols = vfaces_nodes[i].size(); 
        int n_vrtx = n_cols - col00;  
        vector<createpoint> primitive2D(n_vrtx); 
        for(int j=col00, k=0; j< n_cols; j++)
        {
          int vtx = vfaces_nodes[i][j] - 1;
          primitive2D[k++].init( vcoords[vtx], 1);
        }
        gl_model.push_back( primitive2D ); 
      } 
    }
    
    n_vertices = -1; 
    
    cout<<"|_[get_data] Get/GLmodel: ";
    cout<< vfaces_nodes.size() <<"/"; 
    cout<< gl_model.size() <<" "; 
    cout<<"\n";
}


void loader_alya::set_ple_data(const char* cname, int Id)
{
    string fdata = string(cname); 

    // Reading 
    read_log_file Data[3]; 
    Data[0].set_name(fdata+"_coords.alya"); 
    Data[0].run();
    vector< vector<double> > vcoords( Data[0].get_vdata() );
    Data[0].end(); 

    Data[1].set_name(fdata+"_faces_nodes.alya"); 
    Data[1].run();
    vector< vector<double> > vfaces_nodes( Data[1].get_vdata() );
    Data[1].end(); 

    Data[2].set_name(fdata+"_faces_ids.alya"); 
    Data[2].run();
    vector< vector<double> > vfaces_ids( Data[2].get_vdata() );
    Data[2].end(); 

    // Choosing local mesh 
    int col00 = 1; 
    n_elements = 0;
    
    vector<int>  vertex_num; 
    for(int i=0; i<vfaces_ids.size(); i++)
    {
      if(vfaces_ids[i][col00] == Id)
      {
        for(int j=col00; j< vfaces_nodes[i].size(); j++) vertex_num.push_back( vfaces_nodes[i][j] - 0); // base CERO(C++) or ONE(fortran)??
        vertex_num.push_back(-1); // format used by ple?? 
        n_elements++; 
      } 
    }


    vector<double>  vertex_coords;
    for(int i=0; i<vcoords.size(); i++) for(int j=col00; j<vcoords[i].size(); j++) vertex_coords.push_back( vcoords[i][j] ); 

    
    vector< vector<double> > pts; 
    for(int i=0; i<vcoords.size(); i++) 
    {
      vector<double> aux;  
      for(int j=col00; j<vcoords[i].size(); j++) aux.push_back( vcoords[i][j] ); 
      pts.push_back( aux ); 
    }
    
    n_vertices = vcoords.size();

    cout<<"|_[get_data] Got/Chosed/Vrtx: ";
    cout<< vfaces_nodes.size() <<"/"; 
    cout<< n_elements <<"/"; 
    cout<< n_vertices <<" "; 
    cout<<"\n";

    // vector -> ptr
    vertex_num_ptr = new int[ vertex_num.size() ]; 
    memcpy(vertex_num_ptr, vertex_num.data(), sizeof(int)*vertex_num.size() ); 
    vertex_num.clear(); 

    vertex_coords_ptr = new double[ vertex_coords.size() ]; 
    memcpy(vertex_coords_ptr, vertex_coords.data(), sizeof(double)*vertex_coords.size() ); 
    vertex_coords.clear(); 
    
    /*
    vfaces_nodes.clear(); 
    vfaces_ids.clear(); 
    vcoords.clear(); 
    */
}


void loader_alya::print_data()
{
/*
    for(int i=0; i<gl_model.size(); i++)
    {
      for(int j=0; j<gl_model[i].size(); j++)
      {
        cout<< gl_model[i][j].x <<" ";  
        cout<< gl_model[i][j].y <<" ";  
        cout<< gl_model[i][j].z <<" ";  
        cout<<"\n";
      }
      cout<<"\n";
    }

    cout<<"|_[print_data]";
    cout<<"\n";
*/
}


void loader_alya::get_n_vertices()
{
}

void loader_alya::get_n_elements()
{
}

void loader_alya::get_vertex_num_ptr()
{
}

void loader_alya::get_vertex_coords_ptr()
{
}


void loader_alya::draw()
{
}


void loader_alya::draw_quads()
{
}

//=======================================================================||===//
//=======================================================================||===//
