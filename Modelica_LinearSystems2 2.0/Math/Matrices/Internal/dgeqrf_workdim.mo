within Modelica_LinearSystems2.Math.Matrices.Internal;
function dgeqrf_workdim
  "Calculate the optimal size of the WORK array in dgeqrf"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,:];

  output Integer lwork;
  output Integer info;

protected
  Real work[:];

algorithm
  (,,info,work) := LAPACK.dgeqrf(A, -1);
  lwork := integer(work[1]);

end dgeqrf_workdim;
