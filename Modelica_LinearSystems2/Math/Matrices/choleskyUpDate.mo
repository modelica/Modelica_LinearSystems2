within Modelica_LinearSystems2.Math.Matrices;
function choleskyUpDate
  "Compute the cholesky factor Lu according to Au=Lu'*Lu=A + v*v' with A=L'*L"
  extends Modelica.Icons.Function;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real L[:,size(L, 1)] "Cholesky factor";
  input Real v[size(L,1)] "Real vector A+v*v'";
  input Boolean upper=false
    "True if the upper triangle of A is provided and the modified upper triangle will be returned";

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
     i := i+1;
   end while;
 else
    Lud[1,1] := L[1,1]+v[1];
  end if;

  assert(info==0,"Cholesky update failed in choleskyUpDate since updating would not result in a positive definite matrix");

  if upper then //return upper triangle
    for i in 2:n loop
      for j in 1:i-1 loop
        Lud[j,i] := Lud[i,j];
        Lud[i,j] := 0.0;
      end for;
    end for;
  end if;
annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Lud = Matrices.<strong>choleskyUpDate</strong>(L, v);
Lud = Matrices.<strong>choleskyUpDate</strong>(L, v, true);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the rank-1-updated Cholesky factorization <strong>Lud</strong>, with
</p>
<blockquote>
  <strong>Aud</strong> = <strong>Lud</strong>*<strong>Lud</strong><sup>T</sup> =
  <strong>A</strong> +  <strong>v</strong>*<strong>v</strong><sup>T</sup> =
  <strong>L</strong>*<strong>L</strong><sup>T</sup> +  <strong>v</strong>*<strong>v</strong><sup>T</sup>
</blockquote>
<p>
from the input <strong>L</strong>, i.e. the left (lower) Cholesky factor of the original matrix <strong>A</strong>.<br>
The approach is a transformation <strong>H</strong>*[<strong>v</strong>, <strong>L</strong>]' = [<strong>0</strong>, <strong>Lud</strong>]' with orthonormal matrix <strong>H</strong> such, that
</p>
<blockquote>
  [<strong>0</strong>, <strong>Lud</strong>] * [<strong>0</strong>, <strong>Lud</strong>]<sup>T</sup> =
  [<strong>v</strong>, <strong>L</strong>]*<strong>H</strong><sup>T</sup> *<strong>H</strong>*[<strong>v</strong>, <strong>L</strong>]<sup>T</sup> =
  [<strong>v</strong>, <strong>L</strong>]*[<strong>v</strong>, <strong>L</strong>]<sup>T</sup> = <strong>v</strong>*<strong>v</strong><sup>T</sup> + <strong>A</strong>
</blockquote>
<p>
and matrix <strong>Lud</strong> is lower (upper) triangular. The transformation is performed
by n (order of <strong>A</strong>) Givens rotations.
The following sequence illustrates the principle of stepwise transformation of
matrix [v, L]'. Symbol \"*\" represents arbitrary elements. For each step the changed
elements are bold.
</p>
<blockquote><pre>
| v' |    | * * * * |       | 0 <strong>*</strong> <strong>*</strong> <strong>*</strong> |       | 0 0 <strong>*</strong> <strong>*</strong> |       | 0 0 0 <strong>*</strong> |       | 0 0 0 0 |
|    |    | * * * * |       | <strong>*</strong> <strong>*</strong> <strong>*</strong> <strong>*</strong> |       | * * * * |       | * * * * |       | * * * * |
| L' | =  | 0 * * * |  ->   | 0 * * * |  ->   | 0 <strong>*</strong> <strong>*</strong> <strong>*</strong> |  ->   | 0 * * * |  ->   | 0 * * * |
|    |    | 0 0 * * |       | 0 0 * * |       | 0 0 * * |       | 0 0 <strong>*</strong> <strong>*</strong> |       | 0 0 * * |
|    |    | 0 0 0 * |       | 0 0 0 * |       | 0 0 0 * |       | 0 0 0 * |       | 0 0 0 <strong>*</strong> |

</pre></blockquote>
<p>
With the boolean input \"upper\" the user specifies whether the matrix <strong>L</strong>
is lower or upper triangular matrix (left or right Cholesky factor).
If \"upper==true\", the output <strong>Lud</strong> is also upper triangular. Default is \"upper==false\".
</p>
</html>", revisions="<html>
<table border=\"1\" cellspacing=\"0\" cellpadding=\"2\">
  <tr>
    <th>Date</th>
    <th>Author</th>
    <th>Comment</th>
  </tr>
  <tr>
    <td valign=\"top\">2010-05-31</td>
    <td valign=\"top\">Marcus Baur, DLR-RM</td>
    <td valign=\"top\">Realization</td>
  </tr>
</table>
</html>"));
end choleskyUpDate;
