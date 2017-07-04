within Modelica_LinearSystems2.Internal.Streams;
function readMatrixOnFileSize "OBSOLETE- - use Modelica.Utilities.Streams.readMatrixSize instead: Read size of matrix matrixName from file"
  input String fileName "File name";
  input String matrixName "Matrix name";
  output Integer dim[2] "Size of matrix matrixName";
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
}", Icon(graphics={Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={238,46,47},
          lineThickness=0.5,
          pattern=LinePattern.Dash)}));
end readMatrixOnFileSize;
