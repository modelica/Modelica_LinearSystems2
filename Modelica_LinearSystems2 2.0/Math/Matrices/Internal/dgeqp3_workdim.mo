within Modelica_LinearSystems2.Math.Matrices.Internal;
function dgeqp3_workdim
  "Calculate the optimal size of the WORK array in dgeqp3"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,:];

  output Integer lwork;
  output Integer info;

protected
  Real work[:];

algorithm
  (,,,info,work) := LAPACK.dgeqp3(A, -1);
  lwork := integer(work[1]);

end dgeqp3_workdim;
