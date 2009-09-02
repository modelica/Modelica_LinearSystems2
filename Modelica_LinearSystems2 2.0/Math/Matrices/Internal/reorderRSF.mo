within Modelica_LinearSystems2.Math.Matrices.Internal;
function reorderRSF
  "Reorders a real Schur factorization according to a given pattern of the eigenvalues"

  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Boolean iscontinuous;
  input Real T[:,:];
  input Real Q[:,size(T, 2)];
  input Real alphaReal[size(T, 1)]
    "Real part of eigenvalue=alphaReal+i*alphaImag";
  input Real alphaImag[size(T, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

  output Real To[size(T, 1),size(T, 2)];
  output Real Qo[size(T, 1),size(T, 2)];
  output Real wr[size(T, 2)];
  output Real wi[size(T, 2)];

protected
  Integer n=size(T, 2);
  Boolean select[:]=fill(false, size(T, 2));
  Integer i;
algorithm
  if iscontinuous then
    for i in 1:n loop
      if alphaReal[i] < 0 then
        select[i] := true;
      end if;
    end for;
  else
    for i in 1:n loop
      if alphaReal[i]^2 + alphaImag[i]^2 < 1 then
        select[i] := true;
      end if;
    end for;
  end if;

  (To,Qo,wr,wi) := LAPACK.dtrsen(
      "E",
      "V",
      select,
      T,
      Q);

end reorderRSF;
