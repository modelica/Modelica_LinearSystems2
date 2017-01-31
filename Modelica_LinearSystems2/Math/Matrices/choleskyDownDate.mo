within Modelica_LinearSystems2.Math.Matrices;
function choleskyDownDate
  "Compute the cholesky factor Ld according to Ad=Ld'*Ld=A - v*v' with A=L'*L"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real L[:,size(L, 1)] "Cholesky factor";
  input Real v[size(L, 1)] "Real vector A+v*v'";
  input Boolean upper=false "True if the upper triangle of A is provided";

  output Real Ldd[size(L, 1),size(L, 2)]=if upper then symmetric(L) else L
    "Updated Cholesky factor";

protected
  Integer n=size(L, 1);
  Boolean trans=upper;
  Real cvec[size(L, 1)];
  Real svec[size(L, 1)];
  Integer ldL=max(1, n);
  Real a[size(L, 1)]=LAPACK.dtrsv(L, v, upper, trans, false);
  Real q=1 - a*a;
  Integer i;
  Integer info=0;

algorithm
  if q >= 0 then
    if n > 1 then
      q:=sqrt(q);
      for i in n:-1:1 loop
        (cvec[i],svec[i],q) := LAPACK.drotg(q, a[i]);
        if q < 0 then
          q := -q;
          cvec[i] := -cvec[i];
          svec[i] := -svec[i];
        end if;
      end for;

      a:=fill(0,n);
      i := n;
      while i>=1 and info==0 loop
        if Ldd[i, i] < 1e-16 then
          info := -1;
        end if;
        (a[i:n], Ldd[i:n, i]) := LAPACK.drot(a[i:n], Ldd[i:n, i], cvec[i], svec[i]);
        if Ldd[i, i] < 0 then
          Ldd[i:n, i] := -Ldd[i:n, i];
        end if;
        if abs(Ldd[i, i]) < 1e-16 then
          info := -2;
        end if;
        i := i-1;
      end while;
    else
      Ldd[1, 1] := sqrt(L[1, 1]*L[1, 1] - v[1]*v[1]);
    end if;
  else
    info := -3;
//    Modelica.Utilities.Streams.print("q = "+String(q));
  end if;

  assert(info==0,"Cholesky downdate failed in choleskyDownDate since downdating would not result in a positive definite matrix. info = "+String(info));

  if upper then
    for i in 2:n loop
      for j in 1:i-1 loop
        Ldd[j,i] := Ldd[i,j];
        Ldd[i,j] := 0.0;
      end for;
    end for;
  end if;

 annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Ldd = Matrices.<b>choleskyDownDate</b>(L, v);
Ldd = Matrices.<b>choleskyDownDate</b>(L, v, true);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the rank-1-downdated
Cholesky factorization <b>Ldd</b>, with
</p>
<blockquote>
  <b>Add</b> = <b>Ldd</b>*<b>Ldd</b><sup>T</sup> =
  <b>A</b> -  <b>v</b>*<b>v</b><sup>T</sup> =
  <b>L</b>*<b>L</b><sup>T</sup> -  <b>v</b>*<b>v</b><sup>T</sup>
</blockquote>
<p>
from the input <b>L</b>, i.e. the left (lower) Cholesky factor of the
original matrix <b>A</b>. The algortihm is taken from [1].
</p>
<p>
Matrix <b>Ldd</b> is calculated by
</p>
<blockquote>
  [<b>v</b>, <b>Ldd</b>]<sup>T</sup> = <b>H</b> *[<b>0</b>, <b>L</b>]<sup>T</sup>
</blockquote>
<p>
with orthogonal Matrix <b>H</b> such that
</p>
<blockquote>
  <b>v</b>*<b>v</b><sup>T</sup> + <b>Ldd</b>*<b>Ldd</b><sup>T</sup> =
  [<b>v</b>, <b>Ldd</b>] * [<b>v</b>, <b>Ldd</b>]<sup>T</sup> =
  [<b>0</b>, <b>L</b>]*<b>H</b><sup>T</sup> *<b>H</b>*[<b>0</b>, <b>L</b>]<sup>T</sup> =
  [<b>0</b>, <b>L</b>]*[<b>0</b>, <b>L</b>]<sup>T</sup> = <b>L</b>*<b>L</b><sup>T</sup> = <b>A</b>,
</blockquote>
<p>
i.e., by orthogonal transformation
</p>
<blockquote>
  <b>H</b> = <b>H</b>_1*...*<b>H</b>_n.
</blockquote>
<p>
The matrices <b>H</b>_i are Givens matrices computed such that
</p>
<blockquote>
  <b>H</b>_1*<b>H</b>_2*...*<b>H</b>_n*[z, <b>a</b><sup>T</sup>]<sup>T</sup> = [1, 0, ..., 0]<sup>T</sup>,
</blockquote>
<p>
with <b>a</b> is the solution of
</p>
<blockquote>
  <b>L</b>*<b>a</b> = <b>v</b>
</blockquote>
<p>
and
</p>
<blockquote>
  z = ||<b>a</b>||.
</blockquote>
<p>
The following sequence illustrate the principle of calculating the <b>H</b>_i, starting with <b>H</b>_n
</p>
<blockquote><pre>
|z|       |z|       |z|       |z|
|a|  H_3  |a|  H_2  |a|  H_1  |0|
|a| ----> |a| --->  |0| --->  |0|
|a|       |0|       |0|       |0|
</pre></blockquote>
<p>
Note, that the z and a are different in each column.
It is shown in [1] that this algorithm results in the modified Cholesky factor <b>Ldd</b>.
</p>
<p>
With the boolean input \"upper\" the user specifies whether the matrix <b>L</b> is lower
or upper triangular matrix (left or right Cholesky factor).
If \"upper==true\", the output <b>Ldd</b> is also upper triangular. Default is \"upper==false\".
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Dongarra J. J., Bunch J. R., Moler G. B., Stewart G.W. (1987):</dt>
<dd> <b>LINPACK Users' Guide</b>.
     Society for Industrial Mathematics.<br>&nbsp;</dd>
</dl>
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
end choleskyDownDate;
