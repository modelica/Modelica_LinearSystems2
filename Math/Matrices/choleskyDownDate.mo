within Modelica_LinearSystems2.Math.Matrices;
function choleskyDownDate
  "Compute the cholesky factor Ld according to Ad=Ld'*Ld=A - v*v' with A=L'*L"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real L[:,size(L, 1)] "Cholesky factor";
  input Real v[size(L, 1)] "Real vector A+v*v'";
  input Boolean upper=false "True if the upper triangle of A is provided";

  output Real Lud[size(L, 1),size(L, 2)]=if upper then symmetric(L) else L
    "Updated Cholesky factor";

protected
  Integer n=size(L, 1);
  Boolean trans=upper;
  Real cvec[size(L, 1)];
  Real svec[size(L, 1)];
  Integer ldL=max(1, n);
  Real vv[size(L, 1)]=LAPACK.dtrsv(L, v, upper, trans, false);
  Real q=1 - vv*vv;
  Real lii;
  Integer i;
  Integer info=0;

algorithm
  if q > 0 then
    if n > 1 then
      q:=sqrt(q);
      for i in n:-1:1 loop
        lii := Lud[i, i];
        (cvec[i],svec[i],q) := LAPACK.drotg(q, vv[i]);
        if q < 0 then
          q := -q;
          cvec[i] := -cvec[i];
          svec[i] := -svec[i];
        end if;
      end for;

      vv:=fill(0,n);
      i := n;
      while i>=1 and info==0 loop
        if Lud[i, i] < 1e-16 then
//     Modelica.Utilities.Streams.print("ii<0");
          info := -1;
        end if;
        (vv[i:n], Lud[i:n, i]) := LAPACK.drot(vv[i:n], Lud[i:n, i], cvec[i], svec[i]);
        if Lud[i, i] < 0 then
          Lud[i:n, i] := -Lud[i:n, i];
        end if;
        if abs(Lud[i, i]) < 1e-16 then
//    Modelica.Utilities.Streams.print("ii=0");

          info := -1;
        end if;
        i := i-1;
      end while;
    else
      Lud[1, 1] := sqrt(L[1, 1]*L[1, 1] - v[1]*v[1]);
    end if;
  else
//  Modelica.Utilities.Streams.print("q<0");
    info := -1;
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
end choleskyDownDate;
