within Modelica_LinearSystems2.Math.Matrices;
function flipud "Flip the columns of a matrix in up/down direction"
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,:] "Matrix to be flipped";
  output Real Aflip[size(A, 1),size(A, 2)]=A[{i for i in size(A,1):-1:1},:]
    "Flipped matrix";
algorithm
end flipud;
