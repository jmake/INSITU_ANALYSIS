#ifndef OCYL_FILTER_H 
#define OCYL_FILTER_H 

#include "vtktools.hpp" // Analysis, PWriter1  
#include "getzeros.hpp" // GetZeros 

#include <vtkRegularPolygonSource.h>
#include <vtkAlgorithmOutput.h> 
#include <vtkAlgorithm.h>
#include <vtkProbeFilter.h>

#include <vtkArrayCalculator.h> 
#include <vtkGradientFilter.h> 

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
class AddGradients
{
  public :
 ~AddGradients(){};

  AddGradients(vtkAlgorithm* target)
  {
    Target = target->GetOutputPort();  
  };

  vtkArrayCalculator* SetTarget(
                                std::string result, 
                                std::string function, 
                                std::vector<std::string> Scalars,
                                std::vector<std::string> Vectors = std::vector<std::string>()  
                               ) 
  {
    vtkArrayCalculator *calculator = vtkArrayCalculator::New(); 
    calculator->SetInputConnection( Target ); // target->GetOutputPort() );

    for(int i=0; i<Scalars.size(); i++) calculator->AddScalarArrayName( Scalars[i].c_str() ); 
    for(int i=0; i<Vectors.size(); i++) calculator->AddVectorArrayName( Vectors[i].c_str() ); 
    calculator->SetResultArrayName( result.c_str() );
    calculator->SetFunction( function.c_str() ); 
    calculator->Update();

    Target = calculator->GetOutputPort();
    return calculator;  
  }  

  vtkGradientFilter* GetVorticity(std::string key, std::string name="Vorticity")
  {
    vtkGradientFilter *gradients = vtkGradientFilter::New(); 
    gradients->SetInputArrayToProcess(0, 0, 0, 0, key.c_str());
    gradients->SetInputConnection( Target );
    gradients->ComputeGradientOff();
    gradients->ComputeVorticityOn(); 
    gradients->SetVorticityArrayName( name.c_str() );
    gradients->Update(); 
 
    Target = gradients->GetOutputPort();
    return gradients;  
  }

  vtkAlgorithmOutput* GetVtkAlgorithmOutput()  
  {
    return Target; 
  }

  protected :
  vtkAlgorithmOutput *Target ;
  vtkArrayCalculator *Calculator ; 
  vtkGradientFilter  *Gradients; 
};


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
vtkAlgorithm* GetCircle(double Radius, int resolution )
{
  vtkRegularPolygonSource *Target = vtkRegularPolygonSource::New();
  Target->GeneratePolygonOff();
  Target->SetNumberOfSides( resolution );
  Target->SetRadius( Radius );
  Target->SetCenter(0.0, 0.0, 0.0);
  Target->Update();

  return Target;
}


bool CheckRadius(vtkDataObject *dataObject, double r)
{
  assert( dataObject->IsA("vtkUnstructuredGrid") ); //dataObject->PrintSelf(std::cout,vtkIndent(2)); //exit(0);
  //PrintDataObjectName( dataObject );  

  Analysis analysis = Analysis(); // vtktools.hpp   

  if( dataObject->IsA("vtkUnstructuredGrid") ) 
    analysis.SetPoints( vtkUnstructuredGrid::SafeDownCast(dataObject)->GetPoints() ); 
  analysis.PrintArrays(false);
  //analysis.SaveArray("x.y", "Coords");

  std::vector<double> X( analysis.GetComponent("Coords",0) );
  std::vector<double> Y( analysis.GetComponent("Coords",1) );
  std::vector<double> Z( analysis.GetComponent("Coords",2) );

  int rows = X.size(); 
  std::vector<double> mR(rows,0.0);
  for(int i=0; i<rows; i++) mR[i] = sqrt( X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i] );
//for(int i=0; i<rows; i++) std::cout<< i <<" "<< X[i] <<" " << Y[i] <<" "<< Z[i] <<" \n"; 

  std::vector<double>::iterator result = std::min_element(mR.begin(), mR.end());
  int id = std::distance(mR.begin(), result);
  std::cout<<"\t[CheckRadius] Min element at:" << mR[id] <<"= rx"<< mR[id] / r <<" (r="<< r <<") \n"; 

  //exit(0); 
  return true; //mR[id] == r; 
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
class CylinderSurfaceAnalysis : public Analysis // vtktools.hpp 
{
  public :
  CylinderSurfaceAnalysis() : Analysis()
  {
  }


  void FindZeros(std::string xname, std::string yname, int idim=0)
  {
    assert( this->FindArrayName(xname) );  
    assert( this->FindArrayName(yname) );

    std::vector<double> X( Arrays[xname] ); 
    std::vector<double> Y; 

    if(idim) Y = this->GetComponent(yname,idim);
    else     Y = Arrays[yname]; 
    Roots = GetZeros(X,Y); 
    std::cout<<"\tnRoots:"<< Roots.size() <<" \n";

    this->xname = xname;
    this->yname = yname; 
  } 


  void SaveZeros(std::string fname = "roots.dat")
  {
#ifdef USE_BOOST 
    std::vector<double> X( Arrays[xname] ); 
    std::vector<double> Y( Arrays[yname] );  
    boost::math::barycentric_rational<double> b(X.data(), Y.data(), X.size());
    //boost::math::barycentric_rational<double> b(X.data(), Y.data(), X.size());

    std::ofstream myfile;
    myfile.open( fname.c_str() );
    for(int i=0;i<Roots.size();i++) myfile<< Roots[i] <<" "<< b(Roots[i]) <<" \n";
    myfile.close(); 
#endif 
  }


  void Sort()
  {
    Id = this->GetSortId();

    for(It=Arrays.begin(); It!=Arrays.end(); It++)
    {
      std::string name(It->first);
      std::vector<double> Unsorted( It->second );

      int dims = (Unsorted.size()==rows)?(1):( Unsorted.size()/rows );
      std::cout<<"\t'"<< name <<"' ("<< Unsorted.size() <<"/"<< dims <<") ... ";

      std::vector<double*> ptr(dims,NULL);
      for(int i=0; i<dims; i++) ptr[i] = Unsorted.data() + rows * i;

      std::vector<double> Sorted(Unsorted.size(), std::numeric_limits<double>::max());
      for(int i=0; i<rows; i++)
        for(int j=0; j<dims; j++) Sorted[i + j*rows] = ptr[j][ Id[i] ];

      Arrays[name].clear();
      Arrays[name] = Sorted;

      std::cout<< "Sorted! \n";
    }
  }

  std::vector<size_t> GetSortId()
  {
    Coords = Arrays["Coords"];
    double *X = Coords.data() + rows * 0;
    double *Y = Coords.data() + rows * 1;
    double *Z = Coords.data() + rows * 2;

    std::vector<double> Phi(rows,0.0);
    for(int i=0; i<rows; i++) Phi[i] = atan2(Y[i],X[i]);
    Arrays["Phi"] = Phi;
    ArrayNames.push_back("Phi");

    return sort_indexes<double>(Phi);
  }

  protected :
  std::vector<double> Coords;
  std::vector<size_t> Id;
  std::vector<double> Roots; 

  std::string xname;
  std::string yname;
};


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
class FindZeros 
{
  public :
 ~FindZeros(){};

  FindZeros(vtkAlgorithm *source, vtkAlgorithm* target, double radius, int itime=0)
  {
    vtkArrayCalculator *calculator = NULL;

    AddGradients addGradients = AddGradients(source); 
    addGradients.SetTarget(  
                           "VELOC", 
                           "VELOX*iHat + VELOY*jHat + 0*kHat",
                           {"VELOX","VELOY"}  
                          ); 
    addGradients.GetVorticity("VELOC","VORTY"); 

    calculator = 
    addGradients.SetTarget(
                           "mVORTY",
                           "mag(VORTY)", 
                           std::vector< std::string >(),  
                           {"VORTY"}
                          );
    //calculator->PrintSelf(std::cout,vtkIndent(2)); 

    Source   = calculator->GetOutputDataObject(0); 
    assert( CheckRadius(Source, radius) ); //exit(0);  

    Target   = target->GetOutputPort();   
    Probe    = SetInputs(Source, Target);  
    PolyData = GetVtkDataObject(Probe); 

    // TEMPE, VELOX, VELOY, VELOZ, PRESS(NO!), vtkValidPointMask  
    analysis = CylinderSurfaceAnalysis();
    analysis.SetPoints( PolyData->GetPoints() );
    analysis.SetArray(  PolyData->GetPointData(), "VORTY"); // to be used in FindZeros 
    analysis.SetArray(  PolyData->GetPointData(), "vtkValidPointMask"); // interpolation ok? 
    analysis.Sort(); //analysis.PrintArrays(); 
    analysis.FindZeros("Phi", "VORTY", 2);  
  };

  vtkProbeFilter* SetInputs(vtkDataObject *Source, vtkAlgorithmOutput* Target)
  {
    vtkProbeFilter *probe  = vtkProbeFilter::New();
    probe->SetSourceData( Source ); //Source->GetOutputDataObject(0) );
    probe->SetInputConnection( Target ); //Target->GetOutputPort() );
    probe->Update(); 
    //probe->PrintSelf(std::cout,vtkIndent(2));
    return probe; 
  } 

  vtkPolyData* GetVtkDataObject(vtkProbeFilter *probe)
  {
    vtkDataObject *dataObject = probe->GetOutputDataObject(0); //PrintDataObjectName( vtkDataObject );   
    vtkPolyData   *polyData   = vtkPolyData::SafeDownCast(dataObject);  
    return polyData;  
  }

  void SaveVTK(std::string name)
  {
    PWriter1( Probe->GetOutputPort(0), name, Probe->GetOutputDataObject(0)->GetClassName() ); 
  }

  void SavePRO(std::string name, std::string key)  
  {
    std::string fname = name + "." + key;   
    analysis.SaveArray(fname, key);
  }

  void SaveZRS(std::string name)
  {
    std::string fname = name + ".zeros";
    analysis.SaveZeros(fname);
  }

  protected : 
  vtkAlgorithmOutput* Target ;
  vtkProbeFilter *Probe;
  vtkDataObject  *Source; 
  vtkPolyData    *PolyData;
 
  CylinderSurfaceAnalysis analysis;
}; 

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


 
//----------------------------------------------------------------------------
#endif
// J. MIGUEL ZAVALA AKE. 2019AUG18. STOCKHOLM, SWEDEN.
