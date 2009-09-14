within Modelica_LinearSystems2.Internal.Streams;
function readMatrixInternal
  input String fileName "file name";
  input String matrixName "matrix name";
  input Integer m "number of rows";
  input Integer n "numbre of coloumns";
  output Real A[m,n] "matrix to read";
external "C" Modelica_LinearSystem2_readMatrixInternal(
    fileName,
    matrixName,
    A,
    m,
    n);
  annotation (Include="
 
#include <matrixop.h> 
#include <ModelicaUtilities.h>
 
void Modelica_LinearSystem2_readMatrixInternal(const char *filename, const char *matrixname, double* matrix, int rows, int cols) 
{
  int j;
  int i;
  RealArray m;
  m=readMatrix(filename,matrixname,rows,cols); 
  for (j=0;j<rows;j++)
  {
    for (i=0;i<cols;i++)
    {
      matrix[i+j*cols]=m.data[i+j*cols];
    }
  }   
 
//ModelicaFormatMessage(filename);
  return;
}");
end readMatrixInternal;
