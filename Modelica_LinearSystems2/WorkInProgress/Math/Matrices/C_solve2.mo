within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function C_solve2
  "computes the solution to a complex system of linear equations A*X=B, using LU decomposition with partial pivoting and row interchanges"
  import Modelica_LinearSystems2.Math.Complex;
  import Modelica_LinearSystems2.Math.Matrices;

  input Complex A[:,size(A, 1)];
  input Complex B[size(A, 1),:];
  output Complex X[size(A, 1),size(B, 2)];

protected
  Integer info;
algorithm
  (X,info) := Modelica_LinearSystems2.WorkInProgress.Math.LAPACK.zgesv(
                                                                 A, B);
  assert(info == 0, "Solving a linear system of equations with function
\"Matrices.C_solve2\" is not possible, because the system has either
no or infinitely many solutions (A is singular).");
end C_solve2;
