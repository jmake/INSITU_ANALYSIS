#include <string>
#include <vector>
using namespace std;

#ifndef REAL_FILE_H
#define REAL_FILE_H
class read_log_file
{
  public:
    read_log_file(); 
   ~read_log_file();   
    int  set_name(string); 
    void __load_las(); 
    void __find_ranges(); 

    void run(); 
    void end(); 
    void __split_lines(); 
    void __find_logs(); 
    void print_logs(); 

    int get_cols();
    int get_rows(); 
    void get_data_size(int& OUTPUT);
    void get_data(double* OUTPUT); 

    double get_data_ij(int, int);
    vector< vector<double> > get_vdata(); 

  private:
    int     ran; 
    int     n_lines; 
    int     cols;
    int     rows;  
    string* lines; 
    char*   file_name; 
    string *logs_names_array; 
    
    vector< string >          names;
    vector< vector<string> >  split_data;
    vector< vector<double> >  logs;
}; 

#endif 
