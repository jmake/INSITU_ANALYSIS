#ifndef VTK_TOOLS_BASICS_H  
#define VTK_TOOLS_BASICS_H 

#include <iostream> // cin, cout, endl, cerr
#include <vector>   
#include <map>

#include <unistd.h>     // getopt  
#include <assert.h>   
#include <typeinfo>     //for 'typeid' to work  
#include <algorithm>    // std::min_element, std::max_element, std::sort, std::find  

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
#include <vtkXMLPDataReader.h> 
#include <vtkXMLPUnstructuredGridReader.h>
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

#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLWriter.h>

#include <vtkOutlineFilter.h>
#include <vtkDistributedDataFilter.h>

#include <vtkMPI.h>
#include <vtkMPIController.h>
#include <vtkMultiProcessController.h>
#include <vtkMPICommunicator.h>

#include <vtkAppendFilter.h>
#include <vtkVersion.h>
#include <vtkFloatArray.h>

#include <vtkHexahedron.h>
#include <vtkVoxel.h>  
#include <vtkTetra.h>
#include <vtkPixel.h> 
#include <vtkWedge.h>
#include <vtkPyramid.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>

#include <vtkCompleteArrays.h>   

#include <vtkBox.h>
#include <vtkExtractGeometry.h>
#include <vtkImplicitFunction.h>
#include <vtkPlane.h>
#include <vtkCutter.h>

#include <vtkFileOutputWindow.h>

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
void PrintVtkVersion(int print)  
{
  if(print) std::cout<<"\t[PrintVtkVersion] '"<< vtkVersion::GetVTKSourceVersion()  <<"' \n"; 
}


//--------------------------------------------------------------------------||--//
void Warning( std::string fname )
{
  vtkFileOutputWindow *outwin = vtkFileOutputWindow::New();
  outwin->SetFileName( fname.c_str() );    
  outwin->SetInstance(outwin);
}


//--------------------------------------------------------------------------||--//
void
PWriter2(vtkDataObject *Obj, std::string name, vtkMPIController *contr)
{
  assert(!Obj->IsA("vtkMultiBlockDataSet"));

  // 1.0.
  vtkXMLPDataWriter *writer = NULL;
  if(Obj->IsA("vtkUnstructuredGrid" )) writer = vtkXMLPUnstructuredGridWriter::New();
  if(Obj->IsA("vtkPolyData"         )) writer = vtkXMLPPolyDataWriter::New();
  assert(writer);

  // 1.1. 
  writer->SetController(contr);
  writer->SetNumberOfPieces( contr->GetNumberOfProcesses() ); // (nranks);
  writer->SetStartPiece( contr->GetLocalProcessId() ); // (rank);
  writer->SetEndPiece( contr->GetLocalProcessId() ); // (rank);

  // 1.2. 
  if(1) writer->SetDataModeToBinary();
  else  writer->SetDataModeToAscii();

  // 1.3.
  if(1) //  Could not find PPoints element with 1 array 
  {
    vtkNew<vtkCompleteArrays> completeArrays;
    completeArrays->SetController( contr );
    completeArrays->SetInputData( Obj );
    writer->SetInputConnection( completeArrays->GetOutputPort() );
  }
  else
  {
    writer->SetInputData( Obj );
  }

  std::ostringstream oss;
  oss<< name <<"."<< writer->GetDefaultFileExtension();

  writer->SetUseSubdirectory(true);
  writer->SetFileName( oss.str().c_str() );
  writer->Update();

  if(!contr->GetLocalProcessId()) std::cout<<"\t [PWriter2] '"<< name <<"' saved!! \n";
}


//--------------------------------------------------------------------------||--//
class CreateVtkMPIController 
{
  public :
 ~CreateVtkMPIController()
  {
  }


  CreateVtkMPIController(int print=0)
  {
    assert(!contr);

    contr = vtkMPIController::New(); 
    contr->Initialize(NULL, NULL, 1); 
    contr->Barrier();

    nranks = contr->GetNumberOfProcesses();
    rank   = contr->GetLocalProcessId();
    if(!rank && print) std::cout<<"\t[CreateVtkMPIController] "<< nranks <<" Processes Initialized!!\n";  
  }


  CreateVtkMPIController(MPI_Comm lcomm, int print=0)
  {
    assert(!contr);  

    vtkMPICommunicatorOpaqueComm *vtkcomm = new vtkMPICommunicatorOpaqueComm(&lcomm);

    vtkMPICommunicator *comm = vtkMPICommunicator::New();
    comm->InitializeExternal( vtkcomm );

    contr = vtkMPIController::New();
    contr->SetCommunicator( comm );
    contr->Barrier();

    nranks = contr->GetNumberOfProcesses();
    rank   = contr->GetLocalProcessId();
    if(!rank && print) std::cout<<"\t[CreateVtkMPIController] "<< nranks <<" Processes Initialized!!\n";
  }


  void Finalize() // APPLY MPI_Finalize() !!  
  {
    contr->Finalize();  
  }


  vtkMPIController* GetController() 
  {
    assert(contr); 
    return contr; 
  }


  private :
  vtkMPIController *contr = NULL; 
  int rank   = -1;  
  int nranks = -1; 
}; 


//--------------------------------------------------------------------------||--//
template<class VTKARRAY>
class BuildArrays 
{
  public :
 ~BuildArrays()
  {
    this->RemoveArrays();
  }


  BuildArrays() 
  {
    this->RemoveArrays();
  }


  void RemoveArrays()
  {
    for(typename std::map<std::string,VTKARRAY*>::iterator it=this->Arrays.begin(); it!=this->Arrays.end(); ++it)
      it->second->Delete();
      //it->second-Reset(); 

    Arrays.clear();
    Keys.clear();
  }


  void CreateArray(std::string name, int ndata, int ncomps, double *data)
  {
    this->Arrays[name] = VTKARRAY::New();
    this->Arrays[name]->SetName( name.c_str() );
    this->Arrays[name]->SetNumberOfComponents(ncomps);
    this->Arrays[name]->SetNumberOfValues(ndata);
    for(int i=0,k=0; i<ndata; i++) for(int j=0; j<ncomps; j++, k++) this->Arrays[name]->SetValue(k,data[k]);
    this->Keys.push_back(name);
  }


  void UseArray(std::string name, int ndata, int ncomps, double *data)
  {
    this->Arrays[name] = VTKARRAY::New();
    this->Arrays[name]->SetName( name.c_str() );
    this->Arrays[name]->SetNumberOfComponents(ncomps);   
    this->Arrays[name]->SetArray(data, static_cast<vtkIdType>(ndata * ncomps), 1);
    this->Keys.push_back(name);
  }


  VTKARRAY* GetArray(std::string name)
  {
    // SEE : UpdateVTKAttributes, CxxVTKPipelineExample/FEAdaptor.cxx 
    this->ArraysIT = this->Arrays.find(name);
    assert(this->ArraysIT!=this->Arrays.end());
    return this->Arrays[name];
  }


  void PrintArrays()
  {
    for(typename std::map<std::string,VTKARRAY*>::iterator it=this->Arrays.begin();
        it!=this->Arrays.end(); ++it)
    { 
      std::cout<<"\t[PrintArrays] "
	       <<" '"<< it->first <<"'"
	       <<" ("<< it->second->GetNumberOfTuples() 
               <<", "<< it->second->GetNumberOfComponents() <<") "  
	       <<" \n";

      unsigned long kibibytes = 0; 
      kibibytes = it->second->GetActualMemorySize(); // 1024 bytes  
      kibibytes = it->second->GetDataTypeSize() * it->second->GetDataSize();  
      std::string type = it->second->GetDataTypeAsString();  

    } 
  }


  std::vector<std::string> Keys; 

  protected :
  std::map<std::string,VTKARRAY*>            Arrays;
  typename std::map<std::string,VTKARRAY*>::iterator  ArraysIT;
};


//--------------------------------------------------------------------------||--//
template<class intT, class doubleT, class VTKGRID>  
class CreateGrid
{
  public :
 ~CreateGrid()
  {
    this->RemoveGrid(); 
  }; 


  CreateGrid(int dim)
  {
  //this->RemoveGrid();

    this->nDim    = dim;
    this->nCells  = 0;
    this->nPoints = 0;

    this->VTKGrid = VTKGRID::New();
    this->Points  = NULL;   

    this->VTKGridCreated = 0;   
  };


  void RemoveGrid()
  {
    if(this->VTKGrid)
    {
      this->VTKGrid->Delete();
      this->VTKGrid = NULL;
    }

    if(this->Points)
    { 
      this->Points->Delete();
      this->Points = NULL;
    }

    this->nDim = -1;
    this->nCells = -1;
    this->nPoints = -1;

    this->VTKGridCreated = -1;   
  }
 

  void SetPoints(intT npoints, doubleT* pts, int dim)
  {
    assert(dim==3);
    assert(Points == NULL); 
    assert(!VTKGridCreated);  

    Points = vtkPoints::New();
    for(int i=0; i<npoints; i++, this->nPoints++) Points->InsertNextPoint( pts[dim*i+0], pts[dim*i+1], pts[dim*i+2]  );
    VTKGrid->SetPoints(Points);
    Points->Delete();   
    Points = NULL; 
  } 


  template<class T>
  void SetCell(intT nvertices, T *vertices, const char type[], int print)
  {
    assert(!VTKGridCreated);

    std::string t(type);

    int find = 0;
    // 2D : VTK_TRIANGLE(=5); VTK_PIXEL (=8); VTK_QUAD (=9) 
    if(this->nDim==2)
    {
      if(t.find("TRIA" )!=std::string::npos) find = this->SetVtkCell<vtkTriangle  ,T>(nvertices, vertices, print);
      if(t.find("VOXEL")!=std::string::npos) find = this->SetVtkCell<vtkPixel     ,T>(nvertices, vertices, print);
      if(t.find("QUAD" )!=std::string::npos) find = this->SetVtkCell<vtkQuad      ,T>(nvertices, vertices, print);
    }
    // 3D : VTK_TETRA (=10); VTK_VOXEL (=11); VTK_HEXAHEDRON (=12); VTK_WEDGE (=13); VTK_PYRAMID (=14)  
    if(this->nDim==3)
    { 
      if(t.find("TETRA")!=std::string::npos) find = this->SetVtkCell<vtkTetra     ,T>(nvertices, vertices, print);
      if(t.find("VOXEL")!=std::string::npos) find = this->SetVtkCell<vtkVoxel     ,T>(nvertices, vertices, print);
      if(t.find("HEXA" )!=std::string::npos) find = this->SetVtkCell<vtkHexahedron,T>(nvertices, vertices, print);
      if(t.find("WEDGE")!=std::string::npos) find = this->SetVtkCell<vtkWedge     ,T>(nvertices, vertices, print);
      if(t.find("PYRAM")!=std::string::npos) find = this->SetVtkCell<vtkPyramid   ,T>(nvertices, vertices, print);
    } 
    assert(find);
  }

 
  template<class VTKCELL, class T>
  int SetVtkCell(intT nvertices, T *vertices, int print)   
  {
    // X.0. 	  	
    VTKCELL  *cell      = VTKCELL::New();    
    vtkIdType cell_type = cell->GetCellType();   
    vtkIdType cell_size = cell->GetNumberOfPoints(); 
    assert(nvertices==cell_size); 

    // X.0. 
    vtkIdList *vtk_list = vtkIdList::New();  
    for(intT i=0; i<nvertices; i++) vtk_list->InsertNextId( vertices[i] );  
    vtkIdType cell_id = VTKGrid->InsertNextCell(cell_type, vtk_list);  
    vtk_list->Delete();  
    cell->Delete(); 

    // X.0.             
    if(print)
    { 
      std::cout<<"\t[SetVtkCell]"
               <<"   id:"<< cell_id 
               <<" type:"<< cell_type
               <<" size:"<< cell_size
               <<" ["; 
               for(intT i=0; i<nvertices; i++) std::cout<<" "<< vertices[i];
      std::cout<<" ] \n";
    }

    this->nCells += 1;

    return 1;  
  }


  void ON()   
  {
    this->VTKGridCreated = 1;
  }


  void ON(std::string fname)
  {
    this->ON();   
    this->Save(fname); 
  }


  int Created()
  {
    return VTKGridCreated; 
  }


  template<class VTKARRAY>
  void AddPointArray(std::string name, VTKARRAY* array)
  {
    assert(VTKGridCreated);

    int hasArray = VTKGrid->GetPointData()->HasArray(    name.c_str() );
    if(hasArray)   VTKGrid->GetPointData()->RemoveArray( name.c_str() );
    VTKGrid->GetPointData()->AddArray( array );
  }


  void Print()
  {
    VTKGrid->PrintSelf(std::cout,vtkIndent(2));
  }

 
  void Save(std::string fname) 
  {
    assert(VTKGridCreated);
 
    CreateVtkMPIController *MPIController = new CreateVtkMPIController();  
    vtkMPIController *contr = MPIController->GetController(); 

    PWriter2(VTKGrid, fname, contr); 
  }


  VTKGRID       *VTKGrid = NULL;
 
  private :
  intT nDim = -1;
  intT nCells = -1;   
  intT nPoints = -1; 

  vtkPoints      *Points = NULL;  

  int VTKGridCreated = -1; 
}; 



//--------------------------------------------------------------------------||--//
template<class intT, class doubleT, class VTKGRID>
class BuildGrid : public CreateGrid<intT, doubleT, VTKGRID> 
{
  public :
 ~BuildGrid()
  {
    PointProperties.RemoveArrays();  
    CellProperties.RemoveArrays(); 
  }


  BuildGrid(int _dim) : CreateGrid<intT, doubleT, VTKGRID>( _dim ) 
  {
    PointProperties = BuildArrays<vtkDoubleArray>(); 
    CellProperties  = BuildArrays<vtkDoubleArray>();
  }


  void SetPointArray(std::string name, int ndata, int ncomps, double *data)
  {
    PointProperties.CreateArray(name, ndata, ncomps, data); 
  }


  void SetCellArray(std::string name, int ndata, int ncomps, double *data)
  { 
    CellProperties.CreateArray(name, ndata, ncomps, data);                            
  }


  void SetPointArrays()
  {
    for(std::vector<std::string>::iterator it=PointProperties.Keys.begin(); it!=PointProperties.Keys.end(); ++it)
      this->template AddPointArray<vtkDoubleArray>( it[0], PointProperties.GetArray(it[0]) ) ;   
  } 


  private : 
  BuildArrays<vtkDoubleArray> PointProperties; 
  BuildArrays<vtkDoubleArray>  CellProperties;
}; 

/* 
//--------------------------------------------------------------------------||--//
template<class intT, class doubleT>
class BuildVTU : public BuildGrid<intT, doubleT, vtkUnstructuredGrid>
{
  public :
  BuildVTU(int _d) : BuildGrid<intT, doubleT, vtkUnstructuredGrid>(_d) {};


  void ApplyFilter(std::string fname)
  {
    Warning("vtk.log"); 

    // 0.0.  
    CreateVtkMPIController *MPIController = new CreateVtkMPIController();
    vtkMPIController *contr = MPIController->GetController();

    // 1.1. 
    vtkDataObject *obj = this->VTKGrid;  
    assert(obj);  
    PWriter2(obj, fname+"_entire", contr);

    // 1.2. 
    vtkDataObject *extracted = ExtractGeometry(obj, {-1.0,3.0, -1.0,1.0, 0.0,1.0});  // [1]
    PWriter2(extracted, fname+"_extracted", contr);

    // 1.3. 
    vtkDataObject *cutter = Cutter(obj, GetFuntionPlane({2.0,0.0,0.01},{0.0,0.0,1.0}) );
    PWriter2(cutter, fname+"_cutter", contr);

    // 
    // [1] bbox = [xmin,xmax, ymin,ymax, zmin,zmax]
    //
  }

};
*/


//--------------------------------------------------------------------------||--//
//--------------------------------------------------------------------------||--//
#endif
/*
  J. MIGUEL ZAVALA AKE. 2019AUG18. STOCKHOLM, SWEDEN.

  2020JUN15 (vtktools.hpp) : 
  /home/bsc21/bsc21704/z2020_1/REPOSITORY/PAAKAT/TARs/2019AUG09_LIGHT/TESTS/PV560/FORTRAN/E02_BASE

  2020JUN25 (vtktools_basics.hpp) :
  /afs/pdc.kth.se/home/m/miguelza/z2020_1/REPOSITORY/ALYAs/ALYA_2020ABR17/alya/Thirdparties/INSITU/PAAKAT/E04_P4EST11

*/ 
