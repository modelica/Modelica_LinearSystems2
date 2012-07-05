within Modelica_LinearSystems2.Math.Matrices;
function rsf2
  "Computes the real Schur form (RSF) of a square matrix but uses lapack.dgees"
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Matrices.Internal;
  import Modelica_LinearSystems2.Math.Matrices.LAPACK;

  input Real A[:,size(A, 1)];

public
  output Real S[size(A, 1),size(A, 2)];
  output Real QZ[size(A, 1),size(A, 2)];
  output Real alphaReal[size(A, 1)]
    "Real part of eigenvalue=alphaReal+i*alphaImag";
  output Real alphaImag[size(A, 1)]
    "Imaginary part of eigenvalue=(alphaReal+i*alphaImag";

protected
  Integer info;

algorithm
  if size(A, 1) > 1 then
    (S, QZ, alphaReal, alphaImag, info) := Modelica_LinearSystems2.Math.Matrices.LAPACK.dgees(A);
     assert(info == 0, "The output info of LAPACK.dgees should be zero, else if\n
     info < 0:  if info = -i, the i-th argument of dgees had an illegal value\n
     info > 0:  if INFO = i, and i is
               <= N: the QR algorithm failed to compute all the
                     eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
                     contain those eigenvalues which have converged; if
                     JOBVS = 'V', VS contains the matrix which reduces A
                     to its partially converged Schur form.\n");
  else
    S := A;
    if size(A, 1) > 0 then
      QZ := [1];
      alphaReal := {1};
      alphaImag := {0};
    else
      QZ := fill(
        1,
        0,
        0);
      alphaReal := fill(1, 0);
      alphaImag := fill(0, 0);
    end if;
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
(T, Z, alphaReal, alphaImag) = Matrices.<b>rsf2</b>(A)
</pre></blockquote>

<h4>Description</h4>
<p>
Function <b>rsf2</b> (real Schur form) calculates the real Schur form af 
a real square matrix <b>A</b>, i.e.
</p>
<blockquote>
  <b>A</b> = <b>Z</b>*<b>T</b>*<b>Z</b><sup>T</sup>
</blockquote>
<p>
with the real nxn matrices <b>T</b> and <b>Z</b>. <b>Z</b> is an orthogonal matrix. 
<b>T</b> is an block upper triangular matrix with 1x1 and 2x2 blocks in the diagonal. 
The 1x1 blocks contains the real eigenvalues of&nbsp;a. The 2x2 blocks are matrices with 
the conjugated complex pairs of eigenvalues, whereas the real parts of the eigenvalues 
are the elements of the diagonal.
</p>
<p>
The calculation is performed stepwise using lapack.dgees, i.e. using the internal mehtods of balacing and scaling of dgees.
</p>
<p>
See also <a href=\"modelica://Modelica_LinearSystems2.Math.Matrices.rsf\">Math.Matrices.rsf</a>.
</p>

<h4>Example</h4>
<blockquote><pre>
  Real A[3,3] = [1, 2, 3; 4, 5, 6; 7, 8, 9];
  Real T[3,3];
  Real Z[3,3];
  Real alphaReal[3];
  Real alphaImag[3];

<b>algorithm</b>
  (T, Z, alphaReal, alphaImag):=Modelica_LinearSystems2.Math.Matrices.rsf2(A);
//   T = [16.12, 4.9,   1.59E-015;
//        0,    -1.12, -1.12E-015;
//        0,     0,    -1.30E-015]
//   Z = [-0.23,  -0.88,   0.41;
//        -0.52,  -0.24,  -0.82;
//        -0.82,   0.4,    0.41]
//alphaReal = {16.12, -1.12, -1.32E-015}
//alphaImag = {0, 0, 0}
</pre></blockquote>
</html> "));
end rsf2;
