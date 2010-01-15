within Modelica_LinearSystems2.WorkInProgress.Math.Complex.Matrices;
function inv
  "Inverse of a comlex matrix (try to avoid, use function solve(..) instead)"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;

  input Complex A[:,size(A, 1)];
  output Complex invA[size(A, 1),size(A, 2)] "Inverse of matrix A";
protected
  Integer info;
  Integer pivots[size(A, 1)] "Pivot vector";
  Complex LU[size(A, 1),size(A, 2)] "LU factors of A";
algorithm
  (LU,pivots,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgetrf(
                                             A);

  assert(info == 0,
                "Calculating an inverse complex matrix with function
\"Matrices.inv\" is not possible, since complex matrix A is singular.");

  (invA,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgetri(
                                        LU, pivots);

  annotation (Documentation(info=
                             "<html>

</html>"));
end inv;
