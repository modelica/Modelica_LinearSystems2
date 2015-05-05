within Modelica_LinearSystems2.WorkInProgress;
encapsulated function polesAndZerosAsRealMatrix
  "Return poles and zeros of an A,B,C,D system as Real matrix"
  import Modelica_LinearSystems2;

  input Real A[:,size(A,1)] "A-matrix of linear state space system";
  input Real B[size(A,1),:] "B-matrix of linear state space system";
  input Real C[:,size(A,1)] "C-matrix of linear state space system";
  input Real D[size(C,1),size(B,2)] "D-matrix of linear state space system";
  input Boolean poles=true
    "= true, to compute the poles (i.e. the eigenvalues) of (A,B,C,D)"
    annotation (choices(checkBox=true));
  input Boolean zeros=true
    "= true, to compute the (invariant) zeros of (A,B,C,D)"
    annotation (choices(checkBox=true));

  output Real Eigenvalues[:,2]
    "Eigenvalues of A; EigenValues[:,1]: Real part, EigenValues[:,2]: Imaginary part";
  output Real InvariantZeros[:,2]
    "Finite, invariant zeros of linear state space system; size(Zeros,1) <= size(A,1); Zeros[:,1]: Real part, Zeros[:,2]: Imaginary part";
algorithm
  // Determine eigen values
  if poles and size(A, 1) > 0 then
    Eigenvalues := Modelica_LinearSystems2.Math.Matrices.eigenValuesAsRealMatrix(A);
  else
    Eigenvalues := fill(0.0,0,2);
  end if;

  if zeros and size(A,1) > 0 and size(B,2) > 0 and size(C,1) > 0 then
    InvariantZeros := Modelica_LinearSystems2.StateSpace.Internal.invariantZerosWithRealMatrix(A,B,C,D);
  else
    InvariantZeros :=fill(0.0,  0, 2);
  end if;
end polesAndZerosAsRealMatrix;
