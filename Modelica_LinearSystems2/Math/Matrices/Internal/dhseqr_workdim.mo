within Modelica_LinearSystems2.Math.Matrices.Internal;
function dhseqr_workdim
  "Calculate the optimal size of the WORK array in dhseqr"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real H[:,:];

  output Integer lwork;
  output Integer info;

protected
  Real work[:];

algorithm
  if min(size(H, 1), size(H, 2)) > 0 then
    (,,info,,,work) := LAPACK.dhseqr(H, -1);
    lwork := integer(work[1]);
  else
    lwork := 1;
  end if;

end dhseqr_workdim;
