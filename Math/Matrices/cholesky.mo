within Modelica_LinearSystems2.Math.Matrices;
function cholesky
  "Compute the Cholesky factorization of a symmetric positive definte matrix"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real A[:,size(A, 1)];
  input Boolean upper=true "True if the upper triangle of A is provided";

  output Real H[size(A, 1),size(A, 2)]
    "Cholesky factor U or L for A = U'*U or A = L*L'";

protected
  Integer n=size(A,1);
  Integer info;

algorithm
  if size(A, 1) > 0 then
    (H, info) := LAPACK.dpotrf(A, upper);
  else
    H := fill(0,0,0);
    info := 0;
  end if;
  if info<0 then
   assert(info==0,"Cholesky factorization failed in function \"Matrices.cholesky\" due to illegal value of input " +String(info)+" for LAPACK routine DPOTRF");
  else
    assert(info==0,"Cholesky factorization failed in function \"Matrices.cholesky\" since matrix A is not positive definite");
  end if;

  if upper then
    for i in 2:n loop
      for j in 1:i - 1 loop
        H[i, j] := 0.0;
      end for;
    end for;
  else
    for i in 1:n - 1 loop
      for j in i + 1:n loop
        H[i, j] := 0.0;
      end for;
    end for;
  end if;
  annotation (Documentation(info="<html>

</html> ", revisions="<html>
<ul>
<li><i>2010/05/31 </i>
       by Marcus Baur, DLR-RM</li>
</ul>
</html>"));
end cholesky;
