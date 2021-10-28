within Modelica_LinearSystems2.Math.Matrices;
function choleskyDownDate2
  "Compute the cholesky factor Ld according to Ad=Ld'*Ld=A - v*v' with A=L'*L"
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;
  input Real L[:,size(L, 1)] "Cholesky factor";
  input Real v[size(L, 1)] "Real vector A+v*v'";
  input Boolean upper=false "True if the upper triangle of A is provided";

  output Real Ldd[size(L, 1),size(L, 2)]=if upper then symmetric(L) else L
    "Updated Cholesky factor";
  output Integer infoOut;

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
  if q > 1e12 then
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
        if Ldd[i, i] < 1e-12 then
          info := -1;
        end if;
        (a[i:n], Ldd[i:n, i]) := LAPACK.drot(a[i:n], Ldd[i:n, i], cvec[i], svec[i]);
        if Ldd[i, i] < 0 then
          Ldd[i:n, i] := -Ldd[i:n, i];
        end if;
        if abs(Ldd[i, i]) < 1e-12 then
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

//  assert(info==0,"Cholesky downdate failed in choleskyDownDate since downdating would not result in a positive definite matrix. info = "+String(info));
  infoOut := info;

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
Ldd = Matrices.<strong>choleskyDownDate</strong>(L, v);
Ldd = Matrices.<strong>choleskyDownDate</strong>(L, v, true);
</pre></blockquote>

<h4>Description</h4>
<p>
This function computes the rank-1-downdated
Cholesky factorization <strong>Ldd</strong>, with
</p>
<blockquote>
  <strong>Add</strong> = <strong>Ldd</strong>*<strong>Ldd</strong><sup>T</sup> =
  <strong>A</strong> -  <strong>v</strong>*<strong>v</strong><sup>T</sup> =
  <strong>L</strong>*<strong>L</strong><sup>T</sup> -  <strong>v</strong>*<strong>v</strong><sup>T</sup>
</blockquote>
<p>
from the input <strong>L</strong>, i.e. the left (lower) Cholesky factor of the
original matrix <strong>A</strong>. The algortihm is taken from [1].
</p>
<p>
Matrix <strong>Ldd</strong> is calculated by
</p>
<blockquote>
  [<strong>v</strong>, <strong>Ldd</strong>]<sup>T</sup> = <strong>H</strong> *[<strong>0</strong>, <strong>L</strong>]<sup>T</sup>
</blockquote>
<p>
with orthogonal Matrix <strong>H</strong> such that
</p>
<blockquote>
  <strong>v</strong>*<strong>v</strong><sup>T</sup> + <strong>Ldd</strong>*<strong>Ldd</strong><sup>T</sup> =
  [<strong>v</strong>, <strong>Ldd</strong>] * [<strong>v</strong>, <strong>Ldd</strong>]<sup>T</sup> =
  [<strong>0</strong>, <strong>L</strong>]*<strong>H</strong><sup>T</sup> *<strong>H</strong>*[<strong>0</strong>, <strong>L</strong>]<sup>T</sup> =
  [<strong>0</strong>, <strong>L</strong>]*[<strong>0</strong>, <strong>L</strong>]<sup>T</sup> = <strong>L</strong>*<strong>L</strong><sup>T</sup> = <strong>A</strong>,
</blockquote>
<p>
i.e., by orthogonal transformation
</p>
<blockquote>
  <strong>H</strong> = <strong>H</strong>_1*...*<strong>H</strong>_n.
</blockquote>
<p>
The matrices <strong>H</strong>_i are Givens matrices computed such that
</p>
<blockquote>
  <strong>H</strong>_1*<strong>H</strong>_2*...*<strong>H</strong>_n*[z, <strong>a</strong><sup>T</sup>]<sup>T</sup> = [1, 0, ..., 0]<sup>T</sup>,
</blockquote>
<p>
with <strong>a</strong> is the solution of
</p>
<blockquote>
  <strong>L</strong>*<strong>a</strong> = <strong>v</strong>
</blockquote>
<p>
and
</p>
<blockquote>
  z = ||<strong>a</strong>||.
</blockquote>
<p>
The following sequence illustrate the principle of calculating the <strong>H</strong>_i, starting with <strong>H</strong>_n
</p>
<blockquote><pre>
|z|       |z|       |z|       |z|
|a|  H_3  |a|  H_2  |a|  H_1  |0|
|a| ----> |a| --->  |0| --->  |0|
|a|       |0|       |0|       |0|
</pre></blockquote>
<p>
Note, that the z and a are different in each column.
It is shown in [1] that this algorithm results in the modified Cholesky factor <strong>Ldd</strong>.
</p>
<p>
With the boolean input \"upper\" the user specifies whether the matrix <strong>L</strong> is lower
or upper triangular matrix (left or right Cholesky factor).
If \"upper==true\", the output <strong>Ldd</strong> is also upper triangular. Default is \"upper==false\".
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Dongarra J. J., Bunch J. R., Moler G. B., Stewart G.W. (1987):</dt>
<dd> <strong>LINPACK Users' Guide</strong>.
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
end choleskyDownDate2;
