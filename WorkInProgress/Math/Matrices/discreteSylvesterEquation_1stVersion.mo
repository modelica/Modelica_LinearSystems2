within Modelica_LinearSystems2.WorkInProgress.Math.Matrices;
function discreteSylvesterEquation_1stVersion
  "Solution of discrete-time Sylvester equation A*X*B + X = C"
  import Modelica.Math.Matrices;

  input Real A[:,size(A, 1)];
  input Real B[:,size(B, 1)];
  input Real C[size(A, 2),size(B, 1)];
  input Real eps=Matrices.norm(A, 1)*10*Modelica.Constants.eps;

protected
  Integer n=size(A, 1);
  Integer m=size(B, 1);
  Real H[n,n] "hessenberg form  of A, i.e. H=U'AU";
  Real U[n,n] "transformation matrix U for H=U'AU";
  Real S[m,m] "rsf form  of B, i.e. S=Z'BZ";
  Real Z[m,m] "transformation matrix Z for S=Z'BZ";

  Real F[n,m];
  Real R22[n,n];
  Real R11[n,n];
  Integer k;

  Real v1[n];
  Real v2[n];
  Real y[2*size(A, 1)];
  Boolean crit;

public
  output Real X[size(A, 2),size(B, 1)]
    "solution of the discrete Sylvester equation";

algorithm
  X := zeros(n, m);
  k := 1;
  if n > 1 and m > 1 then
    (H,U) := Matrices.hessenberg(A);
    (S,Z) := Matrices.Utilities.rsf(B);
    F := transpose(U)*C*Z;
    while k <= m loop
      v2 := F[:, k] - H*X[:, 1:k - 1]*S[1:k - 1, k];
      crit := if k < m then abs(S[k + 1, k]) < eps else false;
      if (k == m or crit) then
        R22 := S[k, k]*H;
        for i in 1:n loop
          R22[i, i] := R22[i, i] + 1.0;
        end for;
        X[:, k] := Matrices.solve(R22, v2);
        k := k + 1;
      else
        v1 := F[:, k + 1] - H*X[:, 1:k - 1]*S[1:k - 1, k + 1];
        R22 := S[k + 1, k + 1]*H;
        R11 := S[k, k]*H;
        for i in 1:n loop
          R11[i, i] := R11[i, i] + 1.0;
          R22[i, i] := R22[i, i] + 1.0;
        end for;
        y := Matrices.solve([R11,S[k + 1, k]*H; S[k, k + 1]*H,R22], cat(1, v2, v1));
        X[:, k] := y[1:n];
        X[:, k + 1] := y[n + 1:2*n];
        k := k + 2;
      end if;
    end while;

// transform X corresponding to the original form
    X := U*X*transpose(Z);

  elseif n == 1 and m == 1 then
    X[1, 1] := C[1, 1]/(A[1, 1]*B[1, 1] + 1);
  else
    X := fill(
      0,
      0,
      0);
  end if;

  annotation (Documentation(info="<html>
 <h4>Syntax</h4>
<blockquote><pre>
         X = Matrices.<b>continuousLyapunovEquation</b>(A, C);
         X = Matrices.<b>continuousLyapunovEquation</b>(A, C, eps);
</pre></blockquote>
<h4>Description</h4>
<p>
This function computes the solution <b>X</b> of the continuous-time Lyapunov equation
</p>
<blockquote><pre>
 <b>A</b>'*<b>X</b>*<b></b> - <b>X</b> = <b>C</b>.
</pre></blockquote>
<p>
using the Schur method for Lyapunov equations proposed by Bartels and Stewart [1].
<p>
In a nutshell, the problem is reduced to the corresponding problem
<blockquote><pre>
 <b>H</b>*<b>Y</b>*<b>R</b>' - <b>Y</b> = <b>D</b>.
</pre></blockquote>
<p>
with <b>H</b>=<b>U</b>'*<b>A'</b>*<b>U</b> is the RSF of <b>A</b>' and <b>D</b>=<b>U</b>'*<b>C</b>*<b>U</b> and <b>Y</b>=<b>U</b>'*<b>X</b>*<b>U</b>
are the appropriate transformations of <b>C</b> and <b>X</b>. This problem is solved sequently by exploiting the block triangular form of <b>H</b>.
Finally the solution of the the original problem is recovered as <b>X</b>=<b>U</b>*<b>Y</b>*<b>U</b>'.
</p>

<h4>References</h4>
<PRE>
  [1] Bartels, R.H. and Stewart G.W.
      Algorithm 432: Solution of the matrix equation AX + XB = C.
      Comm. ACM., Vol. 15, pp. 820-826, 1972.
</PRE>


</p>
<h4>Example</h4>
<blockquote><pre>
  A = [1, 2,  3,  4;
       3, 4,  5, -2;
      -1, 2, -3, -5;
       0, 2,  0,  6];
  C =  [-2, 3, 1, 0;
        -6, 8, 0, 1;
         2, 3, 4, 5;
        0, -2, 0, 0];
  X = lyapunov(A, C);
  results in:
  X = [1.633, -0.761,  0.575, -0.656;
      -1.158,  1.216,  0.047,  0.343;
      -1.066, -0.052, -0.916,  1.61;
      -2.473,  0.717, -0.986,  1.48]
</pre></blockquote>
<h4>See also</h4>
<a href=\"modelica://Modelica.Math.Matrices.continuousLyapunovEquation\">Matrices.continuousLyapunovEquation</a>


</html>"));
end discreteSylvesterEquation_1stVersion;
