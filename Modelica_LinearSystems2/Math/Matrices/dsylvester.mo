within Modelica_LinearSystems2.Math.Matrices;
function dsylvester
  "Solution of discrete-time Sylvester equation A*X*B + sgn*X = C"
  extends Modelica.Icons.Function;

  import MatricesMSL = Modelica.Math.Matrices;

  input Real A[:,size(A, 1)] "Square matrix A in A*X*B + sgn*X = C";
  input Real B[:,size(B, 1)] "Square matrix B in A*X*B + sgn*X = C";
  input Real C[size(A, 2),size(B, 1)]
    "Rectangular matrix C in A*X*B + sgn*X = C";
  input Boolean AisHess=false "True if A has already Hessenberg form";
  input Boolean BTisSchur=false "True if B' has already real Schur form";
  input Integer sgn=1 "Specifies the sign in A*X*B + sgn*X = C";
  input Real eps=MatricesMSL.norm(A, 1)*10*Modelica.Constants.eps
    "Tolerance";

protected
  Integer n=size(A, 1);
  Integer m=size(B, 1);
  Real H[n,n] "Hessenberg form  of A, i.e. H=U'AU";
  Real U[n,n] "Transformation matrix U for H=U'AU";
  Real S[m,m] "RSF form  of B, i.e. S=Z'BZ";
  Real Z[m,m] "Transformation matrix Z for S=Z'BZ";
  Real F[n,m] "Appropriate transformation of the right side C, F=U'*C*Z";

  Real R22[n,n];
  Real R11[n,n];
  Integer k;

  Real w[n];
  Real g[n];
  Real y[2*n];
  Boolean crit;

public
  output Real X[size(A, 2),size(B, 1)]
    "solution of the discrete Sylvester equation A*X*B + sgn*X = C";

algorithm
  assert(sgn==1 or sgn==-1,"Input sgn in function Math.Matrices.discreteLyapunov() must be 1 or -1, however it is "+String(sgn));
  X := zeros(n, m);
  k := m;

  if n > 1 and m > 1 then
    if AisHess then
      H := A;
      U := identity(n);
      if BTisSchur then
        S := B;
        Z := identity(m);
        F := C;
      else
        (S,Z) := MatricesMSL.realSchur(transpose(B));
        S := transpose(S);
        F := C*Z;
      end if;
    else
      (H,U) := MatricesMSL.hessenberg(A);
      if BTisSchur then
        S := B;
        Z := identity(m);
        F := transpose(U)*C;
      else
        (S,Z) := MatricesMSL.realSchur(transpose(B));
        S := transpose(S);
        F := transpose(U)*C*Z;
      end if;
    end if;

    while k >0 loop

      w := F[:, k] - H*X[:, k + 1:m]*S[k +1:m,k];
      crit := if k > 1 then abs(S[k-1, k]) < eps else false;

      if (k == 1 or crit) then //real eigenvalue in Schur form
        R22 := S[k, k]*H;
        for i in 1:n loop
          R22[i, i] := R22[i, i] + sgn;
        end for;
        X[:, k] := MatricesMSL.solve(R22, w); // solve one column in X for one real eigenvalue
        k := k - 1;
      else // pair of complex eigenvalues, i.e. 2x2 Schur bump
        g := F[:, k-1] - H*X[:, k + 1:m]*S[k+1 :m,k-1];
        R22 := S[k, k]*H;
        R11 := S[k-1, k-1]*H;
        for i in 1:n loop
          R11[i, i] := R11[i, i] + sgn;
          R22[i, i] := R22[i, i] + sgn;
        end for;
        y := MatricesMSL.solve([R11,S[k,  k-1]*H; S[k-1, k]*H,R22], cat(1, g, w));// solve two columns in X for one conjugated complex pole pair
        X[:, k-1] := y[1:n];
        X[:, k] := y[n + 1:2*n];
        k := k - 2;
      end if;
    end while;

// transform X corresponding to the original form
    if not (AisHess and BTisSchur) then
      X := if AisHess then X*transpose(Z) else if BTisSchur then U*X else U*X*transpose(Z);
    end if;

  elseif n == 1 and m == 1 then // simple scalar equation
    X[1, 1] := C[1, 1]/(A[1, 1]*B[1, 1] + sgn);
  else
    X := fill(0, 0, 0);
  end if;

  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
X = Matrices.<strong>dsylvester</strong>(A, B, C);
   or
X = Matrices.<strong>dsylvester</strong>(A, B, C, AisHess, BTisSchur, sgn, eps);
</pre></blockquote>

<h4>Description</h4>
<p>
Function <strong>dsylvester</strong> computes the solution <strong>X</strong> of the discrete-time Sylvester equation
</p>
<blockquote>
  <strong>A</strong>*<strong>X</strong>*<strong>B</strong> + sgn*<strong>X</strong> = <strong>C</strong>.
</blockquote>
<p>
where sgn = 1 or sgn = -1. The algorithm applies the Hessenberg-Schur
method proposed by Golub et al [1]. For sgn = -1, the discrete Sylvester
equation is also known as Stein equation:
</p>
<blockquote>
  <strong>A</strong>*<strong>X</strong>*<strong>B</strong> - <strong>X</strong> + <strong>Q</strong> = <strong>0</strong>.
</blockquote>
<p>
In a nutshell, the problem is reduced to the corresponding problem
</p>
<blockquote>
  <strong>H</strong>*<strong>Y</strong>*<strong>S</strong>' + sgn*<strong>Y</strong> = <strong>F</strong>.
</blockquote>
<p>
with <strong>H</strong>=<strong>U</strong>'*<strong>A</strong>*<strong>U</strong> is the Hessenberg form
of <strong>A</strong> and <strong>S</strong>=<strong>V</strong>'*<strong>B</strong>'*<strong>V</strong> is the real Schur
form of <strong>B</strong>', <strong>F</strong>=<strong>U</strong>'*<strong>C</strong>*<strong>V</strong> and
<strong>Y</strong>=<strong>U</strong>*<strong>X</strong>*<strong>V</strong>' are appropriate transformations
of <strong>C</strong> and <strong>X</strong>. This problem is solved sequently by exploiting
the specific forms of <strong>S</strong> and <strong>H</strong>.
Finally, the solution of the the original problem is recovered as
<strong>X</strong>=<strong>U</strong>'*<strong>Y</strong>*<strong>V</strong>.
</p>
<p>
The boolean inputs \"AisHess\" and \"BTisSchur\" indicate to omit one
or both of the transformation to Hessenberg form or Schur form, respectively,
in the case that <strong>A</strong> and/or <strong>B</strong> have already Hessenberg form
or Schur, respectively.
</p>

<h4><a name=\"References\">References</a></h4>
<dl>
<dt>&nbsp;[1] Golub, G.H., Nash, S. and Van Loan, C.F. (1979):</dt>
<dd> <strong>A Hessenberg-Schur method for the problem AX + XB = C</strong>.
     IEEE Transaction on Automatic Control, AC-24, no. 6, pp. 909-913.<br>&nbsp;</dd>
</dl>

<h4>Example</h4>
<blockquote><pre>
A = [1.0,   2.0,   3.0;
     6.0,   7.0,   8.0;
     9.0,   2.0,   3.0];

B = [7.0,   2.0,   3.0;
     2.0,   1.0,   2.0;
     3.0,   4.0,   1.0];

C = [271.0,   135.0,   147.0;
     923.0,   494.0,   482.0;
     578.0,   383.0,   287.0];

X = discreteSylvester(A, B, C);

// X = [2.0,   3.0,   6.0;
//      4.0,   7.0,   1.0;
//      5.0,   3.0,   2.0];
</pre></blockquote>
</html>",
        revisions="<html>
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
end dsylvester;
