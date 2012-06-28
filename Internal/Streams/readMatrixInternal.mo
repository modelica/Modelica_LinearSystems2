within Modelica_LinearSystems2.Internal.Streams;
function readMatrixInternal "Read matrix matrixName[m,n] from file"
  input String fileName "File name";
  input String matrixName "Matrix name";
  input Integer m "Number of rows";
  input Integer n "Number of coloumns";
  output Real A[m,n] "Matrix to read";
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
}", Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>readMatrixInternal(fileName, matrixName, m, n)</pre></blockquote>

<h4>Description</h4>
<p>Read matrix <code>matrixName</code> from file <code>fileName</code>. The matrix saved in file must be of dimension [m,n].</p>
</html>"));
end readMatrixInternal;
