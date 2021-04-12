within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function C_RQ
  "RQ decomposition of a rectangular complex matrix without column pivoting (A = R*Q)"
  import Complex;

  input Complex A[:,:] "Rectangular matrix with size(A,1) >= size(A,2)";
  output Complex R[size(A, 1),size(A, 1)]
    "Rectangular matrix with orthonormal columns such that Q*R=A[:,p]";
  output Complex Q[size(A, 1),size(A, 2)] "Square upper triangular matrix";

protected
  Integer nrow=size(A, 1);
  Integer ncol=size(A, 2);
  Complex RQ[size(A, 1),size(A, 2)];
  Complex tau[:];

algorithm
  assert(nrow <= ncol, "\nInput matrix A[" + String(nrow) + "," + String(ncol)
     +"] has more rows as columns. This is not allowed when calling Modelica.Matrices.C_RQ(A).");

  if min(nrow,ncol) > 0 then
    (RQ, tau) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgerq2(
                                                                     A);

    for i in 1:nrow loop
      for ii in i:nrow loop
        R[i, ii] := RQ[i, ii+ncol-nrow];
      end for;
    end for;
    Q := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zungrq(
                                                             RQ, tau);
  else
    Q := fill(Complex(1), size(A, 1), 0);
    R := fill(Complex(0), 0,  0);
  end if;
end C_RQ;
