#ifndef VTK_TOOLS_H  
#define VTK_TOOLS_H 

#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>

#include <assert.h>   
#include <algorithm>    // std::min_element, std::max_element, std::sort, std::find  

#include <vtkDataArray.h> 
#include <vtkPointData.h>
#include <vtkPoints.h> 

#include <vtkMPIController.h>
#include <vtkAlgorithmOutput.h>

#include <vtkXMLPDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h> 

#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataObject.h>

#include <vtkUnstructuredGrid.h>
#include <vtkFileOutputWindow.h>


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void Warning( std::string fname )
{
  vtkFileOutputWindow *outwin = vtkFileOutputWindow::New();
  outwin->SetFileName( fname.c_str() );
  outwin->SetInstance(outwin);
}

//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void PrintDataObjectName(vtkDataObject *vtkDataObject)
{
  std::cout<<" GetClassName:'"<< vtkDataObject->GetClassName() <<"' "<<std::endl;
}

vtkDataObject*  ExtractBlock( vtkCompositeDataSet *composite ) 
{
  //vtkCompositeDataSet *composite = reader->GetOutput();
  vtkCompositeDataIterator* iter = composite->NewIterator();
  return iter->GetCurrentDataObject(); 
} 

//--------------------------------------------------------------------------||--//
//----------------------------------------|  FROM : E04_BOOST/nek5k01_01.cxx |--//
std::vector<double>
GetCppArray(vtkDataArray *vtk_array, int* rows=NULL, int* cols=NULL)
{
  assert(vtk_array);

  int nDims = vtk_array->GetNumberOfComponents(); //assert(nDims==1); 
  int nRows = vtk_array->GetNumberOfTuples();

  if(cols) cols[0] = nDims;
  if(rows) rows[0] = nRows;

  std::vector<double> cpp_array(nRows * nDims,0.0);

  /*  C++ Order 
  |--- nDims ---| _ _
  a1  b1 ... y1 z1  |  
  a2  b2 ... y2 z2  |  
  ...              nRows  ==>> [a1 b1 ... y1 z1 ... an bn ... yn zn]    
  an  bn ... yn zn _|_         |---------  nDims * nRows ----------|
  */
  if(0)
  for(int i=0,k=0; i<nRows; i++)
    for(int j=0; j<nDims; j++) cpp_array[k++] = vtk_array->GetComponent(i,j);

  /* Fotran order 
  |--- nDims ---| _ _
  a1  b1 ... y1 z1  |  
  a2  b2 ... y2 z2  |  
  ...              nRows  ==>> [a1 a2 ... an ...       z1 z2 ... zn]    
  a1n bn ... yn zn _|_         |---------  nDims * nRows ----------|
  */
  if(1)
  for(int j=0,k=0; j<nDims; j++)
    for(int i=0; i<nRows; i++) cpp_array[k++] = vtk_array->GetComponent(i,j);

  return cpp_array;
}

//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void PWriter1(vtkAlgorithmOutput *GetOutputPort, std::string name, std::string type)
{
  std::string fname;
  vtkXMLPDataWriter *parallel_writer = NULL;

  if( type == "vtkUnstructuredGrid")
  {
    parallel_writer = vtkXMLPUnstructuredGridWriter::New();
    fname = name + ".pvtu";
  }
  if( type == "vtkPolyData")
  {
    parallel_writer = vtkXMLPPolyDataWriter::New();
    fname = name + ".pvtp";
  }


  vtkMPIController *contr = vtkMPIController::New(); 
  //contr->Initialize(); //&argc, &argv, 1);
  contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

  if(!rank)  
  std::cout
  <<"\t'"<< fname <<"' " 
  <<"rank:"<< rank 
  <<"/"<< nranks
  <<"\n";

  // Create the parallel writer and call some functions
//auto parallel_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  //if( GetOutputPort->IsA("vtkUnstructuredGrid") ) parallel_writer = vtkXMLPUnstructuredGridWriter::New(); 
  //if( GetOutputPort->IsA("vtkPolyData") )    parallel_writer = vtkXMLPPolyDataWriter::New();   
  assert(parallel_writer); 
  parallel_writer->SetInputConnection( GetOutputPort );
  parallel_writer->SetController(contr);
  parallel_writer->SetFileName( fname.c_str() );
  parallel_writer->SetNumberOfPieces(nranks);
  parallel_writer->SetStartPiece(rank);
  parallel_writer->SetEndPiece(rank);
  parallel_writer->SetDataModeToBinary();
  parallel_writer->Update();
  parallel_writer->Write();

  contr->Finalize(1);

/*
  SEE : 
    CxxVTKPipelineExample/vtkCPVTKPipeline.cxx  

  // If process 0 doesn't have any points or cells, the writer may
  // have problems in parallel so we use completeArrays to fill in
  // the missing information.
  vtkNew<vtkCompleteArrays> completeArrays;
  completeArrays->SetInputConnection(threshold->GetOutputPort());

  vtkNew<vtkXMLPUnstructuredGridWriter> writer;
  writer->SetInputConnection(completeArrays->GetOutputPort());
  std::ostringstream o;
  o << dataDescription->GetTimeStep();
  std::string name = this->FileName + o.str() + ".pvtu";
  writer->SetFileName(name.c_str());
  writer->Update();
*/

}

//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
class Analysis
{
  public :
  Analysis(){}
 ~Analysis(){}

  void SetPoints(vtkPoints *points)
  {
    vtkDataArray *dataArray = points->GetData();
    Arrays["Coords"] = GetCppArray(dataArray, &rows, &cols);
    ArrayNames.push_back("Coords");
  }

  void SetArray(vtkPointData  *pointData, std::string key)
  {
    GetArrayNames(pointData);
    FindArrayName(key);

    vtkDataArray *dataArray = pointData->GetArray(key.c_str());
    Arrays[key] = GetCppArray(dataArray);
  }

  void GetArrayNames(vtkPointData *pointData)
  {
    if(ArrayNames.size()==1) // 'Coords' already added...
    {
      int nArrays = pointData->GetNumberOfArrays();
      for(int k=0; k<nArrays; k++) ArrayNames.push_back( pointData->GetArrayName(k) );
    }
  }

  bool FindArrayName(std::string key)
  {
    std::vector<std::string>::iterator it = std::find(ArrayNames.begin(), ArrayNames.end(), key);
    bool found = it != ArrayNames.end();
    if(!found)
    {
      std::cout<<"\t'"<< key <<"' Not Found!! \n" << std::endl;
      for(int k=0; k<ArrayNames.size(); k++) std::cout<<"\t  "<< k <<") '"<< ArrayNames[k] <<"' \n";
      exit(1);
    }
    return found;
  }

  void SaveArray(std::string fname, std::string key)
  {
      std::ofstream myfile;
      myfile.open (fname);
      assert( Arrays.find(key) != Arrays.end() );

      std::vector<double> Array( Arrays[key]  );

      int dims = (Array.size()==rows)?(1):( Array.size()/rows );

      std::vector<double*> ptr(dims,NULL);
      for(int i=0; i<dims; i++) ptr[i] = Array.data() + rows * i;

      for(int i=0; i<rows; i++)
      {
        myfile<< i  <<" ";
        for(int j=0; j<dims; j++)
        {
          myfile<< ptr[j][i] <<" ";
        }
        myfile<<"\n";
      }
      myfile<<" \n";

    myfile.close();
  }

  void PrintArrays(bool entire=true)
  {
    for(It=Arrays.begin(); It!=Arrays.end(); It++)
    {
      std::string name(It->first);
      std::vector<double> Array( It->second );

      int dims = (Array.size()==rows)?(1):( Array.size()/rows );
      std::cout<<"\t[PrintArrays] '"<< name << "': " << Array.size() <<" = "<< Array.size()/dims <<" x "<< dims <<"\n";

      std::vector<double*> ptr(dims,NULL);
      for(int i=0; i<dims; i++) ptr[i] = Array.data() + rows * i;

      if(entire) 
      for(int i=0; i<rows; i++)
      {
        std::cout<< i  <<") ";
        for(int j=0; j<dims; j++)
        {
          std::cout<< ptr[j][i] <<" ";
        }
        std::cout<<"\n";
      }
      std::cout<<" \n";

    } // for 
  }


  std::vector<double> GetComponent(std::string key, int icol)  
  {
    assert( FindArrayName(key) ); 

    std::vector<double> Array( Arrays[key] );  
    int ncols = (Array.size()==rows)?(1):( Array.size()/rows );
 /*
    std::cout<< "Array.size:" << Array.size() <<" \n";   
    std::cout<< "rows:" << rows <<" \n"; 
    std::cout<< "ncols:" << ncols <<" \n"; 
    std::cout<< "icol:" << icol <<" \n";
    std::cout<< "key:" <<  key <<" \n";
exit(0);
 */
    assert( icol + 1 <= ncols ); 

    double *begin = Array.data() + rows * (icol+0);
    double *end   = Array.data() + rows * (icol+1);

    std::vector<double> array(rows, std::numeric_limits<double>::max() ); 
    for(int i=0; i<rows; i++) array[i] = begin[i];//  Array.data() + rows * icol + i;
    return array;  
//    return std::vector<double>(begin,end); 
  }


  protected :
  std::vector<std::string> ArrayNames;
  std::map<std::string, std::vector<double> > Arrays;

  std::map<std::string, std::vector<double> >::iterator It;

  int rows, cols;

};
/*
  Analysis analysis = Analysis();
  analysis.SetPoints( polyData->GetPoints() );
  analysis.SetArray(  polyData->GetPointData(), "Pressure");
  //analysis.SetArray(  polyData->GetPointData(), "Velocity");
*/


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void GetVtkDataObjectMaxMin(vtkDataObject *dataObject)
{
  vtkMPIController *contr = vtkMPIController::New();
  contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

  assert( dataObject->IsA("vtkUnstructuredGrid") ); //dataObject->PrintSelf(std::cout,vtkIndent(2)); //exit(0);
  //PrintDataObjectName( dataObject );  

  Analysis analysis = Analysis(); // vtktools.hpp   

  if( dataObject->IsA("vtkUnstructuredGrid") )
    analysis.SetPoints( vtkUnstructuredGrid::SafeDownCast(dataObject)->GetPoints() );
  //analysis.PrintArrays(false);
  //analysis.SaveArray("x.y", "Coords");

  int idMin = -1; 
  int idMax = -1;
  std::vector<double>::iterator result;
  std::vector< std::vector<double> > range; 
  for(int i=0; i<3; i++)
  {
    std::vector<double> array( analysis.GetComponent("Coords",i) );
    result = std::min_element(array.begin(), array.end());
    idMin  = std::distance(array.begin(), result);

    result = std::max_element(array.begin(), array.end());
    idMax  = std::distance(array.begin(), result);

    range.push_back( {array[idMin], array[idMax]} );
  }

  for(int i=0; i<range.size(); i++)
    std::cout
           <<"\t[GetVtkDataObjectMaxMin] "
           << rank <<"."<< nranks <<":"
	   <<" R"<< i <<" = "
           <<"["<< range[i][0]   
           <<","<< range[i][1] <<"] "
           <<" \n";
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
#endif
// J. MIGUEL ZAVALA AKE. 2019AUG18. STOCKHOLM, SWEDEN.
/*
  NOTES :
    vtkAlgorithmOutput <-- GetOutputPort|SetInputConnection, from vtkAlgorithm  
    vtkDataObject      <-- GetOutputDataObject|SetInputData, from vtkAlgorithm 
*/ 
