#include <iostream> // cin, cout, endl, cerr
#include <sstream>  // stringstream 
#include <cstdlib>  // atoX
#include <cstring>  // strcpy
#include <fstream>
#include <algorithm>
#include "read_file.hpp"
#define uint unsigned int 
#define DELIMITADOR01 "~"

read_log_file::read_log_file()
{
  ran  = -1; 
  cols = -666; 
  rows = -666; 
}


read_log_file::~read_log_file()
{
  if (ran==1)
  {
    //delete[] file_name; 
    //delete[] lines; 
    //delete[] logs_names_array;  
    end(); 
  }
  //cout<<"ok!"<<endl;
}


void read_log_file::end()
{
  if (ran==1)
  {
    //delete[] file_name; 
    delete[] lines;
    split_data.clear(); 
    logs.clear(); 
    //delete[] logs_names_array;  
    ran=0; 
  }
  //cout<<"ok!"<<endl;
}


int read_log_file::set_name(string file_name_in) /* file_name_in -> n_lines */
{
  file_name = new char[file_name_in.size()+1];
  strcpy(file_name, file_name_in.c_str()); // file_name_in -> file_name

  ifstream indata_aux(file_name);          // file_name -> indata_aux 
  if(!indata_aux)
  {
    cerr << "\n    ERROR: File \'"<< file_name <<"\' could not be opened!! \n" << endl;
    exit(1);
  }

  string  line;
  n_lines = 0; 
  while ( indata_aux.good() )  
  {
    getline(indata_aux, line);
    n_lines++;                        // n_lines ??
  }
  indata_aux.close();

  --n_lines;   
  return n_lines; 
}


void read_log_file::run()
{
  __load_las(); 
  __find_ranges(); 
  __split_lines(); // 90% del tiempo 
  __find_logs(); 
  
/*
  __find_logs_names();
  
  ran  = 1; 
  cols = logs_names.size(); 
  rows = logs.size(); 

  logs_names_array = new string[cols]; 
  get_logs_name(logs_names_array); 
*/
}


// PRIVATES /////////////////////////////////////////
void read_log_file::__load_las() // n_lines -> *lines
{
  lines = new string[n_lines];
  
  ifstream indata(file_name);    
  for(int i=0; i<n_lines; i++) getline(indata, lines[i]);
  indata.close();
}


void read_log_file::__find_ranges() /* lines -> dict  */
{
  vector<int> firsts;
  vector<int> lasts;    
  vector<int> size;
  
  for(int i=0; i<n_lines; i++) 
  {
    if(lines[i].find(DELIMITADOR01) != string::npos) // If DELIMITADOR is NOT found, the "npos" is returned.
    {
      string line    = lines[i]; 
      string subline = line.substr(0, line.size()-1);
      
      firsts.push_back(i);
      names.push_back( subline );  
    } 
  }
  firsts.push_back(n_lines); 

  for(uint i=0; i<firsts.size()-1; i++) lasts.push_back( firsts[i+1] );
  for(uint i=0; i<firsts.size()-1; i++)  size.push_back( lasts[i]-firsts[i] );   
  firsts.pop_back(); 
}


void read_log_file::__split_lines() //  *lines -> split_data
{  
  for(int i=0; i<n_lines; i++) 
  {
    string         buff;
    vector<string> row; 
    stringstream   ss( lines[i] ); // Provides an interface to manipulate strings 
                                   // as if they were input/output "streams".                          
                                   // A "stream" can be represented as a source/destination of characters.

    while (ss >> buff) row.push_back(buff); 
    split_data.push_back(row);     // split row!!
  }

  if(split_data.size() != (uint)n_lines)
  {
    cerr << "Error: no lines" << endl;
    exit(1);
  }  
}

	
void read_log_file::__find_logs() /* lines -> logs */
{
  for(int i=0, k=0; i<n_lines; i++)
  {
    // there are not "#" or "~"
    if( (lines[i].find("#")==string::npos)&&(lines[i].find("~")==string::npos) ) 
    { 
      vector<string> row_str( split_data[i] ); // split_data[i] -> row_str 
      
      // row_str -> row_double 
      vector<double> row_double;   
      for(uint j=0; j<row_str.size(); j++) row_double.push_back( atof(row_str[j].c_str()) );
      logs.push_back(row_double); // row_double -> logs  
      k++; 
      
      //if(k > it->second.size) cout << " sizes error.\n";
    } 
  } 
}


void read_log_file::get_data_size(int& total_size) 
{
  int k = 0; 
  for(uint i=0; i<logs.size(); i++) 
  {
    vector<double> row( logs[i] );  
    for(uint j=1; j<row.size(); j++) k++; 
  }
  total_size = k; 

  cout<<"--|Total 1d size: "<< total_size <<endl;
}



void read_log_file::get_data(double* data) 
{
  int k = 0; 
  for(uint i=0; i<logs.size(); i++) 
  {
    vector<double> row( logs[i] );  
    for(uint j=1; j<row.size(); j++) data[k++] = row[j]; 
  }
}


vector< vector<double> > read_log_file::get_vdata()
{
  return logs; 
}



void read_log_file::print_logs() 
{
  for(uint i=0; i<logs.size(); i++) 
  {
    vector<double> row( logs[i] );  
    
    cout<< i+1 <<") ";
    for(uint j=0; j<row.size(); j++) cout<< row[j] <<" "; 
    cout<<endl;
  }
  cout<<"+No rows in logs: "<< logs.size() <<endl;
}


int
read_log_file::get_cols()
{
  vector<int> cols; 
  for(uint i=0; i<logs.size(); i++) 
  {
    vector<double> row( logs[i] );  
    cols.push_back( row.size() ); 
  }

  int max = *std::max_element( cols.data(), cols.data()+cols.size() ); 
  int min = *std::min_element( cols.data(), cols.data()+cols.size() ); 
  if(max!=min) cout<<"WARNING: colums size are diferentes "<< max <<" "<< min <<"!!\n\n"; 
  
  return max; 
}


int
read_log_file::get_rows()
{
  return logs.size();
}


double read_log_file::get_data_ij(int row, int col) 
{
  return logs[row][col]; 
}


/*
using namespace std;
int main()
{
  read_log_file P; 
  
  string fname = "../MeshCplng/Mesh01/mesh01_coords.alya";
  int nline = P.set_name(fname); 
  P.run(); 
  //P.print_logs();
  
  int rows=-1, cols=-1, size=-1;
  P.get_data_size(size); 
  cout<<"\n--|Lines: "<< nline <<" "<< size <<" "<< cols <<endl;

  float *data = new float[size]; 
  P.get_data(data);
  for(int i=0; i<size; i++) cout<< data[i] <<" ";
  cout<<"\n";

  return 0; 
}
*/
