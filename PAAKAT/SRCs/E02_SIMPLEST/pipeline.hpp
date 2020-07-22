#ifndef VTKCPVTKPIPELINE_H
#define VTKCPVTKPIPELINE_H

#include <fstream>      // std::ofstream
#include <string>
#include <vtkCPPipeline.h>
#include <vtkAlgorithmOutput.h>

/*
template <class vtkGridType>
void PWriter(vtkGridType* DataObject=NULL, vtkAlgorithmOutput *GetOutputPort=NULL); 
*/

class vtkCPDataDescription;
class vtkCPPythonHelper;

class vtkCPVTKPipeline : public vtkCPPipeline
{
public:
  static vtkCPVTKPipeline* New();
  vtkTypeMacro(vtkCPVTKPipeline, vtkCPPipeline);
  virtual void PrintSelf(ostream& os, vtkIndent indent);

  virtual int  Finalize(); 

  virtual void Initialize(int outputFrequency, std::string& fileName);

  virtual int RequestDataDescription(vtkCPDataDescription* dataDescription);

  virtual int CoProcess(vtkCPDataDescription* dataDescription);

protected:
  vtkCPVTKPipeline();
  virtual ~vtkCPVTKPipeline();

private:
  vtkCPVTKPipeline(const vtkCPVTKPipeline&) = delete;
  void operator=(const vtkCPVTKPipeline&) = delete;

  int OutputFrequency;
  std::string FileName;

  ofstream myfile; 
  int rank; 

};
#endif
