within Modelica_LinearSystems2.Math.Matrices;
function householderSimilarityTransformation
  "Calculate the similarity transformation SAS of matrix A with householder matrix S = I - 2u*u'"

  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real A[:,size(A, 1)];
  input Real u[size(A, 1)] "householder vector";
  output Real SAS[size(A, 1),size(A, 1)];

protected
  Integer na=size(A, 1);
  Real S[:,:]=-2*matrix(u)*transpose(matrix(u))/(Vectors.length(u)*
      Vectors.length(u));
  Integer i;
algorithm
  for i in 1:na loop
                   //S=I-2u*u'
    S[i, i] := 1.0 + S[i, i];
  end for;
  SAS := S*A*S;

end householderSimilarityTransformation;
