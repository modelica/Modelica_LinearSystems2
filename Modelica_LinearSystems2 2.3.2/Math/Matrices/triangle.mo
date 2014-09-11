within Modelica_LinearSystems2.Math.Matrices;
function triangle "Return the upper/lower triangular part of a square matrix"
  import Modelica_LinearSystems2.Math.Matrices;

  input Real A[:,size(A, 1)] "Square matrix";
  input Boolean upper=true "True for upper triangle to return";
  output Real Tri[size(A, 1),size(A, 2)] "Triangular matrix";
protected
  Integer n=size(A, 1);
algorithm
  if upper then
    for i in 1:n loop
      Tri[i, i:n] := A[i, i:n];
    end for;
  else
    for i in 1:n loop
      Tri[i, 1:i] := A[i, 1:i];
    end for;
  end if;

end triangle;
