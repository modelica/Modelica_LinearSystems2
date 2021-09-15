within Modelica_LinearSystems2.Math.Matrices.Internal;
function symMatMul
  "Calculate the upper triangle of A*B*A'+a*C with B and C symmetric"
  extends Modelica.Icons.Function;
  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,:];
  input Real B[size(A, 2),size(A, 2)];
  input Real C[size(A, 1),size(A, 1)];
  input Boolean add=true "Value is true if a==1, false if a==0";
  output Real M[size(A, 1),size(A, 1)];

protected
  Integer a1=size(A, 1);
  Integer a2=size(A, 2);
  Integer i;
  Integer j;

  Real alpha=1.0;
  Real beta=if add then 1 else 0;

  Real Butri[a2,a2]=B;
  Real Cutri[a1,a1]=C;
  Real Ah[a1,a2];

algorithm
  for i in 1:a2 loop
    Butri[i, i] := B[i, i]/2;
  end for;

  if add then
    for i in 1:a1 loop
      Cutri[i, i] := C[i, i]/2;
    end for;
    for i in 2:a1 loop
      for j in 1:i - 1 loop
        Cutri[i, j] := 0.0;
      end for;
    end for;
  end if;

  Ah := LAPACK.dtrmm(Butri, A, alpha, true, true, false, false);
  M := LAPACK.dgemm(Ah, A, Cutri, alpha, beta, false, true);

// M:= Ah*transpose(A)+Cutri;
  for i in 1:a1 loop
    for j in i:a1 loop
      M[i,j] := M[i,j]+M[j,i];
    end for;
  end for;

  annotation (Documentation(revisions="<html>
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
</html>",        info="<html>
This function is used to efficiently calculate the matrix <strong>X</strong> from equation
<blockquote><pre>
           T
  <strong>X</strong> = <strong>A</strong>*<strong>B</strong>*<strong>A</strong> + <strong>C</strong>.

</pre></blockquote>
with <strong>B</strong> and <strong>C</strong> are symmetric matrices. They hold<blockquote><pre>

   <strong>B</strong> = <strong>B</strong>u + <strong>B</strong>l   and    <strong>C</strong> = <strong>C</strong>u + <strong>C</strong>l,

</pre></blockquote>

where <strong>B</strong>u and <strong>C</strong>u with
<blockquote><pre>
         T               T
  <strong>B</strong>u = <strong>B</strong>l   and   <strong>C</strong>u = <strong>C</strong>l

</pre></blockquote>
are upper triangular matrices. Furthermore, the matrices are defined such that

i.e.,
<blockquote><pre>
          | bij/2  for i = j
  bu,ij = |
          | bij   else

</pre></blockquote>
and cu,ij respectively.<br>
Finally, <strong>X</strong> is given by the sum of a upper triangular matrix and its transposes
<blockquote><pre>
                 T                   T         T                 T              T     T        T
  <strong>X</strong> = <strong>A</strong>*(<strong>B</strong>u+<strong>B</strong>l)*<strong>A</strong> + (<strong>C</strong>u+<strong>C</strong>l) =  <strong>A</strong>*<strong>B</strong>u*<strong>A</strong> + <strong>A</strong>*<strong>B</strong>l*<strong>A</strong> + (<strong>C</strong>u+<strong>C</strong>l) = <strong>A</strong>*<strong>B</strong>u*<strong>A</strong> + <strong>C</strong>u + (<strong>A</strong>*<strong>B</strong>u*<strong>A</strong> + <strong>C</strong>u) =  <strong>E</strong> + <strong>E</strong>

</pre></blockquote>

Since, <strong>X</strong> also has to be symmetric, only the upper triangle of <strong>X</strong> is computed by calculatiing the upper triangle of matrix <strong>E</strong> and adding the upper trinagle of <strong>E</strong>'.<br>
The calculation employs the BLAS functions <strong>dtrmm</strong> and <strong>dgemm</strong>.<br><br>
Note, that only the upper trinagle is calculated. The complete solution could be achieved by the command
<blockquote><pre>
<strong>X</strong> := symmetric(<strong>X</strong>)
</pre></blockquote>
</html>"));
end symMatMul;
