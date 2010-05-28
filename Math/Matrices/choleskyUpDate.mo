within Modelica_LinearSystems2.Math.Matrices;
function choleskyUpDate
  "Compute the cholesky factor Lu according to Au=Lu'*Lu=A + v*v' with A=L'*L"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real L[:,size(L, 1)] "Cholesky factor";
  input Real v[size(L,1)] "Real vector A+v*v'";
  input Boolean upper=false "True if the upper triangle of A is provided";

  output Real Lud[size(L, 1),size(L, 2)]=if upper then symmetric(L) else L
    "Updated Cholesky factor";

protected
  Integer n=size(L,1);
  Real cvec[size(L,1)];
  Real svec[size(L,1)];
  Real vv[size(L,1)]=v;
  Real lii;
  Integer info=0;
  Integer i;

algorithm
  if n>1 then
   i := 1;
   while i<=n and info==0 loop
     lii := Lud[i,i];
     (cvec[i],svec[i],lii) := LAPACK.drotg(lii,vv[i]);
     if lii<0 then
       cvec[i] := -cvec[i];
       svec[i] := -svec[i];
     end if;
     if abs(lii)<1e-16 then
       info:=-1;
     end if;
     (Lud[i:n,i],vv[i:n]) := LAPACK.drot(Lud[i:n,i],vv[i:n],cvec[i],svec[i]);
     Modelica_LinearSystems2.Math.Matrices.printMatrix(Lud,6,"Lud");
     i := i+1;
   end while;
 else
    Lud[1,1] := L[1,1]+v[1];
  end if;

  assert(info==0,"Cholesky downdate failed in choleskyDownDate since the downdating would not result in a positive definite matrix");

  if upper then
    for i in 2:n loop
      for j in 1:i-1 loop
        Lud[j,i] := Lud[i,j];
        Lud[i,j] := 0.0;
      end for;
    end for;
  end if;

  annotation (Documentation(info="<html>

</html> "));
end choleskyUpDate;
