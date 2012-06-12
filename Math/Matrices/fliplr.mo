within Modelica_LinearSystems2.Math.Matrices;
function fliplr "Flip the columns of a matrix in left/right direction"
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,:] "Matrix to be fliped";
  output Real Aflip[size(A, 1),size(A, 2)] = A[:,{i for i in size(A,2):-1:1}]
    "fliped matrix";
algorithm

end fliplr;
