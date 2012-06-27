within Modelica_LinearSystems2.Internal.Streams;
function readMatrixOnFileSize
  input String fileName "file name";
  input String matrixName "matrix name";
  output Integer dim[2] "matrix to read";
external "C" readMatrixSizeEx(
    fileName,
    matrixName,
    dim);
  annotation (Include="
 
#include <matrixop.h> 
 
void readMatrixSizeEx(const char *file, const char *matname, int* dim) 
{
  IntegerArray ms;
  ms=readMatrixSize(file,matname); 
  dim[0]=ms.data[0];
  dim[1]=ms.data[1];
  return;
}");
end readMatrixOnFileSize;
