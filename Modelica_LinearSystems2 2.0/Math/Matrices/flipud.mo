within Modelica_LinearSystems2.Math.Matrices;
function flipud "flip the columns of a matrix in up/down direction"
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,:] "Matrix to be fliped";
  output Real Aflip[size(A, 1),size(A, 2)] "fliped matrix";
algorithm
  for i in 1:size(A, 2) loop
    Aflip[i,:] := A[size(A, 2) + 1 - i,:];
  end for;
end flipud;
