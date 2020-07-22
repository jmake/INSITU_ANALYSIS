#ifndef VTK_TOOLS_H  
#define VTK_TOOLS_H 

#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>

#include <unistd.h>     // getopt  
#include <assert.h>   
#include <typeinfo>  //for 'typeid' to work  
#include <algorithm>    // std::min_element, std::max_element, std::sort, std::find  
#include <numeric>      // std::accumulate


#include <vtkDataArray.h> 
#include <vtkPoints.h> 

#include <vtkAlgorithmOutput.h>
#include <vtkAlgorithm.h>

#include <vtkXMLPDataWriter.h>
#include <vtkXMLPPolyDataWriter.h>
#include <vtkXMLPUnstructuredGridWriter.h> 

#include <vtkCompositeDataIterator.h>
#include <vtkCompositeDataSet.h>
#include <vtkDataObject.h>

#include <vtkUnstructuredGrid.h>

#include <vtkXMLMultiBlockDataReader.h>
#include <vtkPOpenFOAMReader.h>
#include <vtkXMLPDataReader.h> 
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkDataSet.h>

#include <vtkCenterOfMass.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataSetAttributes.h>

#include <vtkPolyData.h> 
#include <vtkCellArray.h>
#include <vtkVertex.h>
#include <vtkDoubleArray.h>

#include <vtkPointSet.h>

#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLWriter.h>

#include <vtkOutlineFilter.h>
#include <vtkDistributedDataFilter.h>

#include <vtkMPI.h>
#include <vtkMPIController.h>
#include <vtkMultiProcessController.h>
#include <vtkMPICommunicator.h>

#include <vtkAppendFilter.h>

#include <vtkXMLDataReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkBox.h>
#include <vtkExtractGeometry.h> 
#include <vtkMultiBlockDataSet.h>
#include <vtkDataSetTriangleFilter.h> 
#include <vtkCleanUnstructuredGrid.h>

#include <vtkStructuredGrid.h>

#include <vtkRegularPolygonSource.h>
#include <vtkTransform.h> 
#include <vtkTransformPolyDataFilter.h>

// TO BE TESTED!! 
//#include <vtkPDistributedDataFilter.h> // Distributed data decomposition (D3) 
#include <vtkPPCAStatistics.h>         // Parallel principal component analysis. 
#include <vtkPOutlineFilter.h>         // An outline consists of the twelve edges of the dataset bounding box (vtk/Filters/Parallel) [1] 
#include <vtkProcessIdScalars.h>       // Add 'ProcessId' scalar 
#include <vtkPProbeFilter.h>           // Probe dataset in parallel 
#include <vtkPDataSetWriter.h>
#include <vtkPDataSetReader.h>         // CanReadFile ?? (Start reading the meta-data pvtk file) [2]  


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
template <class T>
std::map<int, std::vector<int> >
GetDicFromVector3( std::vector<T> vector )
{
  std::map<int, std::vector<int> > dic;
  for(int i=0; i<vector.size(); i++) dic[ vector[i] ].push_back(i);
  return dic;
}


template <class T>
std::vector<int> GetKeys3(std::map< int,std::vector<T> > Dic)
{
  std::vector<int> Keys(Dic.size(), std::numeric_limits<int>::max());

  int i=0;
  typename std::map< int,std::vector<T> >::iterator It;
  for(It=Dic.begin(); It!=Dic.end(); ++It, i++) Keys[i] = It->first;
  return Keys;
}


template <class T>
std::vector<int> GetKeys3(std::map<int,int> Dic)
{
  std::vector<int> Keys(Dic.size(), std::numeric_limits<int>::max());

  int i=0;
  typename std::map<int,int>::iterator It;
  for(It=Dic.begin(); It!=Dic.end(); ++It, i++) Keys[i] = It->first;
  return Keys;
}


template <class T>
std::vector<T> GetValues3(std::map< int,std::vector<T> > Dic)
{
  std::vector<T> Values(Dic.size(), std::numeric_limits<int>::max());

  int i=0;
  typename std::map< int,std::vector<T> >::iterator It;
  for(It=Dic.begin(); It!=Dic.end(); ++It, i++) Values[i] = It->second;
  return Values;
}


template <class T>
std::vector<T> GetValues3(std::map<int,int> Dic)
{
  std::vector<T> Values(Dic.size(), std::numeric_limits<int>::max());

  int i=0;
  typename std::map<int,int>::iterator It;
  for(It=Dic.begin(); It!=Dic.end(); ++It, i++) Values[i] = It->second;
  return Values;
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
vtkMPIController* GetVtkMPIControllerFromMPIComm(MPI_Comm lcomm)
{
  vtkMPICommunicatorOpaqueComm *vtkcomm = new vtkMPICommunicatorOpaqueComm(&lcomm);

  vtkMPICommunicator *comm = vtkMPICommunicator::New();
  comm->InitializeExternal( vtkcomm );

  vtkMPIController *contr = vtkMPIController::New();
  contr->SetCommunicator(comm);
/*
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();
  std::cout<<"[GetVtkMPIControllerFromMPIComm] rank.nranks:'"<< rank+1 <<"."<< nranks <<" \n";
*/
  return contr;  
};


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void PrintDataObjectName(vtkDataObject *vtkDataObject)
{
  std::cout<<" GetClassName:'"<< vtkDataObject->GetClassName() <<"' "<<std::endl;
}


//--------------------------------------------------------------------------||--//
vtkDataObject*  ExtractBlock( vtkCompositeDataSet *composite ) 
{
  //vtkCompositeDataSet *composite = reader->GetOutput();
  vtkCompositeDataIterator* iter = composite->NewIterator();
  return iter->GetCurrentDataObject(); 
} 

//--------------------------------------------------------------------------||--//
//----------------------------------------|  FROM : E04_BOOST/nek5k01_01.cxx |--//
std::vector<double>
GetCppArray(vtkDataArray *vtk_array, int* rows=NULL, int* cols=NULL, std::string order="F")
{
  //assert(vtk_array);
  std::vector<double> cpp_array;  

  if(vtk_array) 
  { 
    int nDims = vtk_array->GetNumberOfComponents(); //assert(nDims==1); 
    int nRows = vtk_array->GetNumberOfTuples();

    if(cols) cols[0] = nDims;
    if(rows) rows[0] = nRows;

  //std::vector<double> cpp_array(nRows * nDims,0.0);
    cpp_array = std::vector<double>(nRows * nDims, 0.0);   

    /*  C++ Order 
    |--- nDims ---| _ _
    a1  b1 ... y1 z1  |  
    a2  b2 ... y2 z2  |  
    ...              nRows  ==>> [a1 b1 ... y1 z1 ... an bn ... yn zn]    
    an  bn ... yn zn _|_         |---------  nDims * nRows ----------|
    */
//    if(0)
    if(order == "C")
    for(int i=0,k=0; i<nRows; i++)
      for(int j=0; j<nDims; j++) cpp_array[k++] = vtk_array->GetComponent(i,j);

    /* Fotran order 
    |--- nDims ---| _ _
    a1  b1 ... y1 z1  |  
    a2  b2 ... y2 z2  |  
    ...              nRows  ==>> [a1 a2 ... an ...       z1 z2 ... zn]    
    a1n bn ... yn zn _|_         |---------  nDims * nRows ----------|
    */
//    if(1)
    if(order == "F")
    for(int j=0,k=0; j<nDims; j++)
      for(int i=0; i<nRows; i++) cpp_array[k++] = vtk_array->GetComponent(i,j);
  }

  return cpp_array;
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void PWriter1(vtkAlgorithmOutput *GetOutputPort, std::string name, std::string type)
{
  vtkMPIController *contr = vtkMPIController::New();
  contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

  vtkXMLPDataWriter *parallel_writer = NULL;

  std::string fname;
  if( type == "vtkUnstructuredGrid")
  {
    parallel_writer = vtkXMLPUnstructuredGridWriter::New();
    //fname = name + ".pvtu";
  }
  if( type == "vtkPolyData")
  {
    parallel_writer = vtkXMLPPolyDataWriter::New();
    //fname = name + ".pvtp";
  }

  std::ostringstream oss;
  oss<< name
//  <<"_"<< std::setfill('0') << std::setw(6) << nranks
  <<"."<< parallel_writer->GetDefaultFileExtension();

  fname = oss.str(); 

  if(!rank) std::cout<<"\t [PWriter1] '"<< fname <<"' saved!! \n";

  // Create the parallel writer and call some functions
//auto parallel_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  //if( GetOutputPort->IsA("vtkUnstructuredGrid") ) parallel_writer = vtkXMLPUnstructuredGridWriter::New(); 
  //if( GetOutputPort->IsA("vtkPolyData") )    parallel_writer = vtkXMLPPolyDataWriter::New();   
  assert(parallel_writer);
  //parallel_writer->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), -1);
  //parallel_writer->UpdateTimeStep(-1.0); 
  parallel_writer->SetInputConnection( GetOutputPort );
  parallel_writer->SetController(contr);
  parallel_writer->SetFileName( fname.c_str() );
  parallel_writer->SetNumberOfPieces(nranks);
  parallel_writer->SetStartPiece(rank);
  parallel_writer->SetEndPiece(rank);
  parallel_writer->SetDataModeToBinary();
  parallel_writer->Update();
  parallel_writer->Write();

  //contr->Finalize(1);

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
void PVtkObjectWriter(vtkDataObject *vtkObj, std::string name, std::string type)
{
  vtkMPIController *contr = vtkMPIController::New();
  contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

  vtkXMLPDataWriter *parallel_writer = NULL;

  std::string fname;
  if( type == "vtkUnstructuredGrid")
  {
    parallel_writer = vtkXMLPUnstructuredGridWriter::New();
    //fname = name + ".pvtu";
  }
  if( type == "vtkPolyData")
  {
    parallel_writer = vtkXMLPPolyDataWriter::New();
    //fname = name + ".pvtp";
  }
  assert(parallel_writer);

  std::ostringstream oss;
  oss<< name
//<<"_"<< std::setfill('0') << std::setw(6) << nranks
  <<"."<< parallel_writer->GetDefaultFileExtension();

  fname = oss.str();

  if(!rank) std::cout<<"\t [PWriter1] '"<< fname <<"' saved!! \n";

  // Create the parallel writer and call some functions
//auto parallel_writer = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
  //if( GetOutputPort->IsA("vtkUnstructuredGrid") ) parallel_writer = vtkXMLPUnstructuredGridWriter::New();
  //if( GetOutputPort->IsA("vtkPolyData") )    parallel_writer = vtkXMLPPolyDataWriter::New();
  assert(parallel_writer);
  //parallel_writer->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), -1);
  //parallel_writer->UpdateTimeStep(-1.0);
//parallel_writer->SetInputConnection( GetOutputPort );
  parallel_writer->SetInputData( vtkObj );
  parallel_writer->SetController(contr);
  parallel_writer->SetFileName( fname.c_str() );
  parallel_writer->SetNumberOfPieces(nranks);
  parallel_writer->SetStartPiece(rank);
  parallel_writer->SetEndPiece(rank);

//parallel_writer->SetDataModeToBinary();
  parallel_writer->SetDataModeToAscii(); 

  parallel_writer->Update();
  parallel_writer->Write();

}



//--------------------------------------------------------------------------||--//
void VtkObjectWriter(vtkDataObject *vtkObj, std::string name, std::string type)
{
  vtkXMLWriter *writer = NULL; 

  std::string fname; 
  if(type == "vtkPolyData")         writer = vtkXMLPolyDataWriter::New();
  if(type == "vtkUnstructuredGrid") writer = vtkXMLUnstructuredGridWriter::New(); 
  assert(writer);

  std::ostringstream oss;
  oss<< name
  <<"."<< writer->GetDefaultFileExtension();

  fname = oss.str();
  writer->SetFileName( fname.c_str() );
  writer->SetInputData( vtkObj );
  // 
  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  //writer->SetDataModeToAscii();
  writer->Write(); 

  std::cout<<"\t [VtkObjectWriter] '"<< fname <<"' saved!! \n";
}


//--------------------------------------------------------------------------||--//
vtkDataObject* PReader(std::string name, std::string type, int parallel=0)
{
  vtkMPIController *contr = vtkMPIController::New();
  //contr->Initialize(); //&argc, &argv, 1);
  contr->Initialize(NULL, NULL, 1); // initializedExternally == 1;
  int nranks = contr->GetNumberOfProcesses();
  int rank   = contr->GetLocalProcessId();

//  vtkDataSet *x = NULL; 
  vtkUnstructuredGrid *dataSet = vtkUnstructuredGrid::New(); 
  vtkDataObject *obj = NULL; //vtkDataObject::New(); 
//  vtkDataSet *dataSet = vtkDataSet::New();

  if(!rank)
  { 
    vtkXMLPDataReader *parallel_reader = NULL;

    if(type == "vtu") parallel_reader = vtkXMLPUnstructuredGridReader::New();
    if(type == "vtp") parallel_reader = vtkXMLPPolyDataReader::New(); 
    assert(parallel_reader);

    parallel_reader->SetFileName(name.c_str());
    parallel_reader->Update(); //parallel_reader->Print(std::cout);

    obj = parallel_reader->GetOutputDataObject(0); //obj->Print(std::cout);
    assert(obj);

    if(obj->IsA("vtkUnstructuredGrid")) dataSet = vtkUnstructuredGrid::SafeDownCast(obj);
  //if(obj->IsA("vtkPolyData")        ) dataSet = vtkPolyData::SafeDownCast(obj);
  } 

  assert(dataSet);
  if(!rank) std::cout<<"[PReader] vtuNumberOfPoints:"<< dataSet->GetNumberOfPoints() <<" \n";

  if(parallel)
  {
    vtkDistributedDataFilter *distFilter = vtkDistributedDataFilter::New();  
    distFilter->SetInputData( dataSet );
    distFilter->SetController( contr );
    distFilter->Update();  //distFilter->GetOutputDataObject(0)->PrintSelf(std::cout,vtkIndent(2));

    obj = distFilter->GetOutputDataObject(0); 
    std::cout<<"[PReader"<< rank <<"] vtuNumberOfPoints:"<< vtkUnstructuredGrid::SafeDownCast(obj)->GetNumberOfPoints() <<" \n";
  } 

  return obj;
}


//--------------------------------------------------------------------------||--//
vtkDataObject* Reader(std::string fname, std::string type)
{
  vtkXMLDataReader *r = NULL; 
  if(type == "vtu") r = vtkXMLUnstructuredGridReader::New();   
  if(type == "vts") r = vtkXMLStructuredGridReader::New();   
  if(type == "vtp") r = vtkXMLPolyDataReader::New();  
//  if(type == "vtm") r = vtkXMLMultiBlockDataReader::New();

  if(fname.find("vtm"  ) != std::string::npos) std::cout<<"[Reader] Use VtmReader!! \n";
  if(fname.find("foam" ) != std::string::npos) std::cout<<"[Reader] Use FoamReader!! \n"; 
  assert(r);

  r->SetFileName( fname.c_str() );
  r->Update(); // r->Print(std::cout);

  std::cout<<"\t [Reader] '"<< fname <<"'!! \n";
  return r->GetOutputDataObject(0); 
} 


vtkDataObject* FoamReader(std::string fname)
{
  vtkPOpenFOAMReader *reader = vtkPOpenFOAMReader::New();
  reader->SetFileName( fname.c_str() );
//reader->SetController( contr );
  reader->Update(); // reader->Print(std::cout); 

  std::cout<<"\t [Reader] '"<< fname <<"'!! \n";
  return reader->GetOutputDataObject(0); 
//  vtkDataObject *obj = TetrahedralizeMultiblock( reader->GetOutputDataObject(0) );
//  return vtkUnstructuredGrid::SafeDownCast( obj );
} 


vtkDataObject* VtmReader(std::string fname)
{
  vtkXMLMultiBlockDataReader *reader = vtkXMLMultiBlockDataReader::New();
  reader->SetFileName( fname.c_str() );
  reader->Update(); // reader->Print(std::cout);

  std::cout<<"\t [Reader] '"<< fname <<"'!! \n";
  return reader->GetOutputDataObject(0);
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
vtkDataObject* AddProcessId(vtkDataObject* obj, vtkMultiProcessController *ctr)   
{
  vtkProcessIdScalars *d = vtkProcessIdScalars::New();
  d->SetInputData( obj );

  d->SetController( ctr ); 
  d->SetScalarModeToPointData(); 

  d->Update();

  return d->GetOutputDataObject(0);    
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
vtkDataObject* GetCircle(double Radius, int resolution )
{
  vtkRegularPolygonSource *Target = vtkRegularPolygonSource::New();
  Target->GeneratePolygonOff();
  Target->SetNumberOfSides( resolution );
  Target->SetRadius( Radius );
  Target->SetCenter(0.0, 0.0, 0.0);
  Target->Update();

//  return Target;
  return Target->GetOutputDataObject(0);
}


//--------------------------------------------------------------------------||--//
vtkPolyData *GetOutline(vtkDataObject *vtkObj)
{
  vtkOutlineFilter *outline = vtkOutlineFilter::New();
  outline->SetInputData( vtkObj );
  outline->Update();

  return vtkPolyData::SafeDownCast( outline->GetOutputDataObject(0) );
}


//--------------------------------------------------------------------------||--//
vtkDataObject* GetTransform(vtkDataObject* obj, std::vector<double> scale, std::vector<double> translate)   
{
  vtkTransform *transform = vtkTransform::New(); 
  transform->Scale( &scale[0] );
  transform->Translate( &translate[0] ); 
  //transform->RotateY(0.0); 

  vtkTransformPolyDataFilter *tf = vtkTransformPolyDataFilter::New(); 
  tf->SetInputData( obj ); 
  tf->SetTransform( transform );
  tf->Update(); 

//  vtkDataObject * dummy; 
//  return dummy;  
  return tf->GetOutputDataObject(0);
}



//--------------------------------------------------------------------------||--//
std::vector<double> GetPBBox(vtkDataObject *vtkObj, MPI_Comm lcomm, int debug=0)
{
  assert(lcomm!=MPI_COMM_NULL);

  vtkMPIController *c = GetVtkMPIControllerFromMPIComm(lcomm); 

  vtkPOutlineFilter *outline = vtkPOutlineFilter::New();
  outline->SetInputData( vtkObj );
  outline->SetController( c ); 
  outline->Update();

  vtkPolyData *vtp = vtkPolyData::SafeDownCast( outline->GetOutputDataObject(0) ); 
  double *bounds = vtp->GetBounds();
  std::vector<double> bbox(bounds, bounds+6);

  c->Broadcast(&bbox[0], bbox.size(), 0); 
  
  if( debug )
  {  
    std::cout<<"\t[GetPOutline."<< debug <<"] BBox=[ ";
    for(int i=0; i<bbox.size(); i++) std::cout<< bbox[i] <<" ";
    std::cout<<"] \n";
  } 

  return bbox; 
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
std::vector<std::string>
Argvs(int argc, char* argv[])
{
  std::vector<std::string> files;
  int c ;
  while((c = getopt(argc, argv, "a:b:")) != -1)
  {
    switch(c)
    {
      case 'a':
        if(optarg) files.push_back(optarg); //avalue = optarg;
        break;
      case 'b':
        if(optarg) files.push_back(optarg); //bvalue = optarg;
        break;
    }
  }

  std::cout<<"[Argvs] nFiles:"<<files.size() <<" ";
  for(int i=0; i<files.size(); i++) std::cout<<" '"<< files[i] <<"' ";
  std::cout << '\n' ;

  return files;
}


//--------------------------------------------------------------------------||--//
std::vector<double> 
GetCenterOfMassFilter(vtkDataObject *vtkObj)
{
  vtkDataSet *ptData = NULL;
  if(vtkObj->IsA("vtkUnstructuredGrid")) ptData = vtkUnstructuredGrid::SafeDownCast(vtkObj); //->GetPointData();
  if(vtkObj->IsA("vtkPolyData"        )) ptData = vtkPolyData::SafeDownCast(vtkObj); //;->GetPointData();
  assert(ptData);

  std::vector<double> CoM; //(3,std::numeric_limits<int>::max());  
  if( ptData->GetNumberOfPoints() )
  {	   
    vtkCenterOfMass *centerOfMassFilter = vtkCenterOfMass::New();
    centerOfMassFilter->SetInputData( vtkObj );
    centerOfMassFilter->SetUseScalarsAsWeights(false);
    centerOfMassFilter->Update();

    double *ptr = centerOfMassFilter->GetCenter(); 
    CoM = std::vector<double>(ptr,ptr+3);    
  }

  return CoM; 
}


//--------------------------------------------------------------------------||--// 
std::vector<double> 
GetPointDataArray(vtkDataObject *vtkObj, std::string key)
{
   vtkPointData *data = NULL;
   if(vtkObj->IsA("vtkUnstructuredGrid")) data = vtkUnstructuredGrid::SafeDownCast(vtkObj)->GetPointData();
   if(vtkObj->IsA("vtkPolyData"        )) data = vtkPolyData::SafeDownCast(vtkObj)->GetPointData();
   assert(data);

   std::map<std::string,int> Dic;
   for(int i=0; i<data->GetNumberOfArrays(); i++) Dic[data->GetArrayName(i)] = i; 
   if( Dic.find(key.c_str()) == Dic.end() ) std::cout<<"[GetPointDataArray] WARNING:'"<< key.c_str()  <<"' no found!! \n";

   vtkDataArray *Array = data->GetArray( key.c_str() ); //assert(Array);
   std::vector<double> CppArray;
   if(Array) CppArray = GetCppArray(Array);

   return CppArray;
}


std::vector<double>
GetCellDataArray(vtkDataObject *vtkObj, std::string key)
{ 
   vtkCellData *data = NULL;
   if(vtkObj->IsA("vtkUnstructuredGrid")) data = vtkUnstructuredGrid::SafeDownCast(vtkObj)->GetCellData();
   if(vtkObj->IsA("vtkPolyData"        )) data = vtkPolyData::SafeDownCast(vtkObj)->GetCellData();
   assert(data);

   std::map<std::string,int> Dic;
   for(int i=0; i<data->GetNumberOfArrays(); i++) Dic[data->GetArrayName(i)] = i; 
   if( Dic.find(key.c_str()) == Dic.end() ) std::cout<<"[GetPointDataArray] WARNING:'"<< key.c_str()  <<"' no found!! \n";

   vtkDataArray *Array = data->GetArray( key.c_str() ); //assert(Array);
   std::vector<double> CppArray;
   if(Array) CppArray = GetCppArray(Array);

   return CppArray;
}


//std::vector< std::vector<double> >  
//std::vector<double>
void 
GetPtsArray(vtkDataObject *vtkObj) 
{
  vtkPointSet *data = NULL;
  if(vtkObj->IsA("vtkUnstructuredGrid")) data = vtkUnstructuredGrid::SafeDownCast(vtkObj); //->GetPointData();
  if(vtkObj->IsA("vtkPolyData"        )) data = vtkPolyData::SafeDownCast(vtkObj); //;->GetPointData();
  assert(data);

/*
  std::vector<double> CppArray;
  vtkPoints *points = data->GetPoints();
  if(points)
  {
    vtkDataArray *Array = points->GetData(); //if(!dataArray) { std::cout<<" !dataArray"; exit(0); } 
    if(Array) CppArray = GetCppArray(Array); 
  }
*/
  std::vector<double> CppArray;

  vtkPoints *Pts = data->GetPoints();

  int nPts = Pts->GetNumberOfPoints();   

  std::vector<double> X(nPts);
  std::vector<double> Y(nPts);
  std::vector<double> Z(nPts);
  for(int i=0; i<nPts; i++)
  {
    double *pt = Pts->GetPoint(i);  
    X[i]=pt[0]; Y[i]=pt[1]; Z[i]=pt[2];    
  }

  std::map<std::string, std::vector<double> > Coords;
  Coords["X"] = X; 
  Coords["Y"] = Y;
  Coords["Z"] = Z;

//  return CppArray;
}



std::vector<double>
GetPtsArray(vtkDataObject *vtkObj, int dim)
{
  vtkPointSet *data = NULL;
  if(vtkObj->IsA("vtkUnstructuredGrid")) data = vtkUnstructuredGrid::SafeDownCast(vtkObj); //->GetPointData();
  if(vtkObj->IsA("vtkPolyData"        )) data = vtkPolyData::SafeDownCast(vtkObj); //;->GetPointData();
  assert(data);

  vtkPoints *Pts = data->GetPoints();
  int nPts = Pts->GetNumberOfPoints();

  std::vector<double> CppArray(nPts);  
  for(int i=0; i<nPts; i++) CppArray[i] = Pts->GetPoint(i)[dim]; 

  return CppArray;
}


//--------------------------------------------------------------------------||--//
vtkPolyData *VectorDic2VtkVertices( std::map<int, std::vector<double> > Dic )   
{
  vtkPoints         *Pts = vtkPoints::New(); 
  vtkCellArray    *Cells = vtkCellArray::New(); 
  vtkDoubleArray *Arrays = vtkDoubleArray::New(); 

  Arrays->SetName("Ids"); 
  Arrays->SetNumberOfComponents(1); 
  Arrays->SetNumberOfValues( Dic.size() ); 

  typename std::map< int,std::vector<double> >::iterator It;  
  for(It=Dic.begin(); It!=Dic.end(); ++It) 
  {
    int idx = Pts->InsertNextPoint( It->second.data() ); 
    Arrays->InsertValue(idx,It->first);

    vtkVertex *vertex = vtkVertex::New(); 
    vertex->GetPointIds()->SetId(0,idx); // Is '0' because only one vertex exist? 
    Cells->InsertNextCell(vertex);   
  }

  vtkPolyData *Polydata = vtkPolyData::New();
  Polydata->SetPoints(Pts);   
  Polydata->SetVerts(Cells);   
  Polydata->GetPointData()->SetScalars(Arrays); 
  return Polydata; 
}


//--------------------------------------------------------------------------||--//
void AddArrayFromDic(vtkPolyData *Poly, std::string key, std::map<int,double> Dic)
{
  vtkDoubleArray *Arrays = vtkDoubleArray::New(); 
  Arrays->SetName( key.c_str() ); 
  Arrays->SetNumberOfValues( Dic.size() );
  Arrays->SetNumberOfComponents( 1 );   

  int ii=0;
  typename std::map<int,double>::iterator It; 
  for(It=Dic.begin(); It!=Dic.end(); ++It) Arrays->SetValue(ii++,It->second); 

  Poly->GetPointData()->AddArray(Arrays); 
}


void AddArrayFromDic(vtkPolyData *Poly, std::string key, std::map<int, std::vector<double> > Dic, int size)
{
  vtkDoubleArray *Arrays = vtkDoubleArray::New();
  Arrays->SetName( key.c_str() );
  Arrays->SetNumberOfTuples( Dic.size() );
  Arrays->SetNumberOfComponents( size );

  int ii=0; 
  typename std::map<int, std::vector<double> >::iterator It;
  for(It=Dic.begin(); It!=Dic.end(); ++It)  
  {
    assert( It->second.size()==size );   
    Arrays->InsertTuple(ii++,It->second.data());
  //Arrays->InsertTuple(It->first,It->second.data()); // :(   
  } 

  Poly->GetPointData()->AddArray(Arrays);
}


//--------------------------------------------------------------------------||--//
vtkPolyData *CombinePolyData(std::vector< vtkSmartPointer<vtkDataObject> > Polys)
{
  vtkPolyData* Merged = NULL; //vtkPolyData::New(); 

  if( Polys.size() ) 
  {
    vtkAppendPolyData *append = vtkAppendPolyData::New();
    for(int i=0; i<Polys.size(); i++) 
    { 
      vtkPolyData *p = vtkPolyData::SafeDownCast(Polys[i]);
      if(p)
      {
//        std::cout<<"\t[CombinePolyData] "<< i <<" "<< p->GetPoints()->GetNumberOfPoints() <<" \n";
        append->AddInputData( p );
      }
      else
      {
//        std::cout<<"\t[CombinePolyData] "<< i <<" 0!"; 
      }

    } 
    append->Update();  

    vtkCleanPolyData *cleaned = vtkCleanPolyData::New(); 
    cleaned->SetInputConnection( append->GetOutputPort() );
    cleaned->Update(); 

    Merged = vtkPolyData::SafeDownCast( cleaned->GetOutputDataObject(0) );
    if(Merged) std::cout<<"\t[CombinePolyData] "<< Merged->GetPoints()->GetNumberOfPoints() <<" \n"; 
  }

  return Merged;  
}

 
//--------------------------------------------------------------------------||--//
template <class T>
std::vector<T> VtkAllGatherV(std::vector<T> sendBuffer, MPI_Comm comm)
{
  int nranks = -1; 
  int rank   = -1; 
  MPI_Comm_size(comm, &nranks);
  MPI_Comm_rank(comm, &rank  );
//MPI_Datatype Tmpi = U;  

  int sendLength = sendBuffer.size(); 
  std::vector<int> recvLengths(nranks, -1);
 
  MPI_Allgather(&sendLength,     1, MPI_INTEGER,  
		&recvLengths[0], 1, MPI_INTEGER, comm);
 
  std::vector<int>              offsets(nranks + 1, 0);
  for (int i = 0; i < nranks; ++i) offsets[i + 1] = recvLengths[i] + offsets[i] ;  

  std::vector<T> recvBuffer( offsets[nranks] );
  if( recvBuffer.size() ) 
  {
    MPI_Allgatherv(&sendBuffer[0], sendLength,                   MPI_INTEGER, 
                   &recvBuffer[0], &recvLengths[0], &offsets[0], MPI_INTEGER, comm);  
/*
    std::cout<<"\t[VtkAllGatherV."<< rank <<"] "<< sendLength <<" [";
    for(int i=0; i<recvBuffer.size(); i++) std::cout<<" "<< recvBuffer[i] ;
    std::cout<<" ] \n";
*/ 
  }

  return recvBuffer;
}
 

/* 
//
// For some reason that I dont understand, the c->AllGather blocks the processors!!
// Apparently, MPI_Allgatherv works correctly. 
//
template <class T>
std::vector<T> VtkAllGatherV(std::vector<T> sendBuffer, MPI_Comm comm)
{
  vtkMPIController *c = GetVtkMPIControllerFromMPIComm(comm); 
  int nranks = c->GetNumberOfProcesses();
  int rank   = c->GetLocalProcessId();

  vtkIdType                sendLength = static_cast<vtkIdType>(sendBuffer.size());
  std::vector< vtkIdType > recvLengths(nranks, 0);
  c->AllGather(&sendLength, &recvLengths[0], 1);
 
  std::vector< vtkIdType >         offsets(nranks + 1, 0);
  for (int i = 0; i < nranks; ++i) offsets[i + 1] = recvLengths[i] + offsets[i];

  std::vector<T> recvBuffer( offsets[nranks] );
  if( recvBuffer.size() )
  {
    c->AllGatherV(&sendBuffer[0],   
                  &recvBuffer[0], sendLength, &recvLengths[0], &offsets[0]);
  }  

  c->Finalize();
  return recvBuffer;
}
// */



//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
template <class T>
vtkDataObject* MergeVtkDic(std::map<int,T*> Dic)
{
  vtkAppendFilter *append = vtkAppendFilter::New();
  T *dummy = T::New();

  typename std::map<int,T*>::iterator It;

  for(It=Dic.begin(); It!=Dic.end(); ++It) append->AddInputData( (Dic.size())?(It->second):(dummy) );
  append->Update();
  assert( append->GetOutputDataObject(0)->IsA("vtkUnstructuredGrid") );

  return append->GetOutputDataObject(0);
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
template<class T>
void IterateDic(std::map<int,T> D)
{    
  for(typename std::map<int,T>::iterator It=D.begin(); It!=D.end(); ++It) F(&D); 
};



//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
vtkDataObject* ExtractGeometry(vtkDataObject *Obj, std::vector<double> bbox)
{
  vtkBox *b = vtkBox::New();
  b->SetBounds( bbox.data() );

  vtkExtractGeometry *eg = vtkExtractGeometry::New();
  eg->SetInputData( Obj );
  eg->SetImplicitFunction(b);
  eg->ExtractInsideOn();
  eg->ExtractBoundaryCellsOn();
  eg->Update();
  return eg->GetOutputDataObject(0);
}


vtkDataObject* Tetrahedralize(vtkDataObject *Obj)
{
  vtkDataSetTriangleFilter *getTria = vtkDataSetTriangleFilter::New();
  getTria->SetInputData(Obj);
  getTria->Update();
  return getTria->GetOutputDataObject(0);
}


vtkDataObject* TetrahedralizeMultiblock(vtkDataObject *Obj)
{
  assert( Obj->IsA("vtkMultiBlockDataSet") );

  std::vector<double> bbox{0.21,0.26, -0.105,-0.055, 0.015, 0.065};
  //std::for_each(bbox.begin(), bbox.end(),[](double &x){x*=1.1;}); 

  vtkDataObject *Region = Tetrahedralize( ExtractGeometry(Obj,bbox) ) ;

  vtkMultiBlockDataSet *mb = vtkMultiBlockDataSet::SafeDownCast( Region ); // ( Obj );  
  vtkCompositeDataIterator* iter = mb->NewIterator();
  vtkAppendFilter  *append = vtkAppendFilter::New();
  for(iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextItem())
  {
    vtkDataObject* inputBlock = iter->GetCurrentDataObject();
    append->AddInputData(inputBlock);
  }
  append->Update(); // dirty 

  vtkCleanUnstructuredGrid *clean = vtkCleanUnstructuredGrid::New();
  clean->SetInputData( append->GetOutputDataObject(0) );
  clean->Update();

  vtkDataObject *merged = clean->GetOutputDataObject(0);
  return merged;
}


vtkDataObject* ProbeWith(vtkDataObject *source, vtkDataObject* input)
{
/*
  vtkDataSet *input = NULL;
  if( Obj->IsA("vtkPolyData")         ) input = vtkPolyData::SafeDownCast(Obj);
  if( Obj->IsA("vtkUnstructuredGrid") ) input = vtkUnstructuredGrid::SafeDownCast(Obj);
  assert(input);  
*/ 
  vtkProbeFilter *probeFilter = vtkProbeFilter::New();
  probeFilter->SetSourceData( source );  
  probeFilter->SetInputData(  input  );  
  probeFilter->Update();

  //probeFilter->GetOutputDataObject(0)->Print(std::cout);  
  return probeFilter->GetOutputDataObject(0);
}



//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
std::vector<int> GetVerticesArray(vtkDataObject *Obj, int print=1)
{
  // 
  // Vertices = [ cell1, cell2, ..., cellN ] ; celli = [vertex1, vertex2, ..., vertexM];     
  // Types    = [ type1, type2, ..., typeN ] ; typei = celli type   
  // Sizes    = [ size1, size2, ..., sizeN ] ; sizei = len(celli) 
  // Cumus    = [     0, sum(Sizes[:1]), sum(Sizes[:2]), ..., sum(Sizes[:N])]
  // 
  //                cell1      cell2  cell3         
  // Vertices = [ 4 19 23; 11 7 13 4;   ...    ]
  // Types    = [    tria       quad    ...    ]
  // Sizes    = [       3          4    ...    ]
  // Cumus    = [       0          3      7 ...]
  // 
  vtkDataSet *ds = NULL;
  if( Obj->IsA("vtkPolyData"        ) ) ds = vtkPolyData::SafeDownCast(Obj);
  if( Obj->IsA("vtkStructuredGrid")   ) ds = vtkStructuredGrid::SafeDownCast(Obj);
  if( Obj->IsA("vtkUnstructuredGrid") ) ds = vtkUnstructuredGrid::SafeDownCast(Obj);
  assert(ds);

  int    nCells = ds->GetNumberOfCells();
  int nVertices = 0;  
  std::vector<int> Types(nCells,-1);
  std::vector<int> Sizes(nCells,-1);
  for(vtkIdType id=0; id<nCells; id++)
  {
    vtkIdList *icell = vtkIdList::New();
    ds->GetCellPoints(id,icell);
    nVertices += icell->GetNumberOfIds(); 

    Sizes[id] = icell->GetNumberOfIds();    
    Types[id] =    ds->GetCellType(id);
  }
  int nTypes = GetDicFromVector3<int>(Types).size();

  std::vector<int> Vertices(nVertices,-1);
  for(int id=0,k=0; id<nCells; id++)
  {
    vtkIdList *icell = vtkIdList::New();
    ds->GetCellPoints(id,icell);
    for(int j=0; j<icell->GetNumberOfIds(); j++ ) Vertices[k++] = icell->GetId(j);
  }

  if(print)
  {
    std::vector<int> Cumus(nCells+1,0); 
    std::partial_sum( Sizes.begin(), Sizes.end(), &Cumus[1] );   

    for(int i=0,k=0; i<Types.size(); i++)
    {
      std::cout<<" "<< i
	       <<" "<< Types[i]  
               <<" "<< Sizes[i]   
               <<" "<< Cumus[i]
	       <<" ["; 

      for(int j=Cumus[i]; j<Cumus[i+1]; j++, k++) std::cout<<" "<< Vertices[j] ;  

      std::cout<<" ] \n";
    } 

    Cumus.clear();  
  }

  Types.clear(); 
  Sizes.clear(); 

  std::cout<<"[GetVerticesArray] nCells:"<< nCells 
                             <<" nTypes:"<< nTypes 
			     <<" \n"; 

  assert(nTypes == 1); // Only type=10 (VTK_TETRA) is considered!!  
  return Vertices;
}


std::vector<double> GetPointsArray(vtkDataObject *Obj, int print=1)
{
  // 
  // points  = [x1,y1,z1, x2,y2,z2, ..., xN,yN,zN]  
  //
  vtkPoints *points = NULL;
  if( Obj->IsA("vtkPolyData")         ) points = vtkPolyData::SafeDownCast(Obj)->GetPoints();
  if( Obj->IsA("vtkStructuredGrid")   ) points = vtkStructuredGrid::SafeDownCast(Obj)->GetPoints();  
  if( Obj->IsA("vtkUnstructuredGrid") ) points = vtkUnstructuredGrid::SafeDownCast(Obj)->GetPoints(); 
  assert(points); 

  int nPts=0, nCols=0;
  std::vector<double> array(0);
  if(points)
  {
    vtkDataArray *dataArray = points->GetData(); assert(dataArray);
    array = GetCppArray(dataArray, &nPts, &nCols, "C");
  }

  if(print)
  {
    for(int i=0,k=0; i<nPts; i++)
    {
      std::cout<<" "<< i <<" ["; 
      for(int j=0; j<nCols; j++, k++) std::cout<<" "<< array[k] ; 
      std::cout<<" ] \n"; 
    } 
  } 

//std::cout<<"[GetPointsArray] nPts:"<< nPts <<" \n";
  return array;
}


//--------------------------------------------------------------------------||--//
std::vector<double> GetPleppPoints(vtkDataObject *Obj, int dims, int print)
{
  // points  = [x1,y1,z1, x2,y2,z2, ..., xN,yN,zN] 
  std::vector<double> vtk_points = GetPointsArray(Obj, 0); //print); 

  int vtk_dims = 3;
  int vtk_npts = vtk_points.size() / vtk_dims; // dim(vtk) = 3; 

  std::vector<double> ple_points(vtk_npts * dims); 
  for(int i=0,k=0; i<vtk_points.size(); i+=vtk_dims)
    for(int j=i; j<i+dims; j++) ple_points[k++] = vtk_points[j] ; 

  if(print)
  for(int i=0,k=0; i<vtk_points.size(); i+=vtk_dims)
  {
    std::cout<<" "<< i <<" [";
    for(int j=i; j<i+dims    ; j++) std::cout<<" "<< ple_points[k++] ;
    std::cout<<" ] \n"; 
  }

  vtk_points.clear(); 
  vtk_points = ple_points; 
  ple_points.clear();
// /* 
  std::cout<<"[GetPleppPoints] nPts:"<< vtk_npts   
                             <<" "<< vtk_points.size() / dims  
                             <<" \n";
// */
  assert(vtk_npts * dims == vtk_points.size()); 
  return vtk_points;   
}


std::vector<int> GetPleppCells(vtkDataObject *Obj, int print)
{
  // 1.0. 
  // VTK FORMAT (C++ index [0,N-1] ): 
  //  
  //   connectivities = [cell1, cell2, ..., cellN], celli = [vi1, vi2,..., viM]  
  //                  = [v11,v12,...,v1M, v21,v22,...,v2M, ..., vN1,vN2,...,vNM]
  //                       0   1 ...   M    0   1 ...   M         0   1 ...   M
  // 
  //   EXAMPLE : 
  //          2 
  //          +        VTK_TRIANGLE(=5)
  //         / \       cell = [0,1,2], M=3 
  //        /   \
  //       +-----+     
  //     0         1
  //
  std::vector<int> connectivities = GetVerticesArray(Obj, print);

  // 2.0.
  // PLEPP FORMAT (FORTRAN index [1,N]):
  std::for_each(connectivities.begin(), connectivities.end(), [](int &c){c+=1;}); // Fortran style for ple  
  for(auto c:connectivities) assert(c);  
//for(auto c:connectivities) std::cout<< c <<" "; std::cout << '\n';

  return connectivities;
}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
void Dummy()  
{


}


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
#endif
// J. MIGUEL ZAVALA AKE. 2019AUG18. STOCKHOLM, SWEDEN.
/*
  NOTES :
    vtkAlgorithmOutput <-- GetOutputPort|SetInputConnection, from vtkAlgorithm  
    vtkDataObject      <-- GetOutputDataObject|SetInputData, from vtkAlgorithm 


[1] 
  https://github.com/Kitware/VTK/blob/master/Filters/Parallel/vtkPOutlineFilterInternals.cxx
[2]
  https://github.com/Kitware/VTK/blob/master/Parallel/MPI/Testing/Cxx/TestPProbe.cxx


*/ 
