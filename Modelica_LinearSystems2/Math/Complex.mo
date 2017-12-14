within Modelica_LinearSystems2.Math;
operator record Complex "Record defining a Complex number"

  Real re "Real part of complex number" annotation(Dialog);
  Real im "Imaginary part of complex number" annotation(Dialog);

  encapsulated package Examples
    "Package of examples to demonstrate the usage of complex numbers"

    import Modelica;
    extends Modelica.Icons.ExamplesPackage;

    function addTwoComplexNumbers "Show how to add two complex numbers"
      import Complex;
      import Modelica.ComplexMath.j;
      import Modelica.Utilities.Streams;

      output Boolean ok;
    protected
      Complex c1=2+3*j;
      Complex c2=3+4*j;
      Complex c3;
    algorithm

      Streams.print("c1 = " + String(c1));
      Streams.print("c2 = " + String(c2));

      c3 := c1 + c2;
      Streams.print("c3 = c1 + c2 = " + String(c3));
      ok := true;
      annotation (Documentation(info="<html>
<p>In this example there are defined two complex numbers and the sum of these numbers is given on output.</p>
</html>"));
    end addTwoComplexNumbers;
  end Examples;

  encapsulated package Vectors
    "Package of functions operating on vectors of complex numbers"
    extends Modelica.Icons.Package;
    import Modelica;

    function print "Print vector"
      extends Modelica.Icons.Function;

      import Modelica;
      import Modelica.Utilities.Streams.print;
      import Complex;
      input String name="" "Name of complex vector";
      input Complex c[:] "Complex vector to be printed";

      input String fileName=""
        "Print to fileName; empty fileName prints to the log window";

    algorithm
      print(name + " =", fileName);
      for i in 1:size(c, 1) loop
        print("   " + String(c[i]), fileName);
      end for;
    end print;

    function printHTML
      "Print complex vector as HTML in sorted form (vector is expected to have pure real and/or conjugate complex values (so poles and/or zeros)"
      extends Modelica.Icons.Function;

      import Modelica;
      import Modelica.Utilities.Files;
      import Modelica.Utilities.Streams.readFile;
      import Modelica.Utilities.Streams.print;
      import Complex;
      import Modelica_LinearSystems2;

      input Complex c[:] "Complex vector to be printed";
      input String heading="Zeros" "Heading above the table";
      input String name="zeros" "Heading of value column";
      input Boolean sort=true
        "= true, if values are sorted, otherwise no sorting";
      input Boolean ascending = true
        "= true if ascending order, otherwise descending order for sorting";
      input Boolean sortFrequency=true
        "= true, if sorting is first for imaginary then for real value, otherwise sorting is for absolute value";
    protected
      Integer nc = size(c,1);
      Complex cSorted[size(c,1)];
      Integer cIndex[size(c,1)];
      String tempFile = "TemporaryForPrint.html";
      Integer nReal;

      encapsulated function printTable
        "Print the table with eigenvalues in html format on file"
        import Modelica;
        import Modelica.Utilities.Strings;
        import Modelica.Utilities.Streams.print;
        import Modelica_LinearSystems2;
        import Modelica_LinearSystems2.Internal.Eigenvalue;
        import Complex;

        input Complex systemZeros[:];
        input Integer nReal;
        input String tempFile;
      protected
        Integer nz=size(systemZeros, 1);

        String number;
        Real timeConstant;
        Real freq;
        Real damp;

      algorithm
        for i in 1:nReal loop
          // Build eigenvalue number

          number := String(
                  i,
                  minimumLength=7,
                  leftJustified=false);
          timeConstant := if abs(systemZeros[i].re) > 10*Modelica.Constants.eps
             then 1/abs(systemZeros[i].re) else 1/(10*Modelica.Constants.eps);

          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
             + number + " </td>" + "\n  <td> &nbsp; " + String(systemZeros[i].re,
            format="14.4e") + " </td>" + "\n  <td> &nbsp; " + String(
            timeConstant, format="9.4f") + " </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>\n</tr>",
            tempFile);

        end for;

        for i in nReal + 1:2:nz loop
          number := String(i) + "/" + String(i + 1);
          number := Strings.repeat(max(0, 7 - Strings.length(number))) + number;

          // Determine frequency and number of corresponding zero
          (freq,damp) := Modelica_LinearSystems2.Math.Complex.frequency(systemZeros[i]);

          print(
            "<tr style=\"background-color:white\">\n  <td style=\"text-align:left\"> &nbsp; "
             + number + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(systemZeros[i].re, format="14.4e") + " &plusmn; " +
            String(systemZeros[i].im, format="12.4e") + "j </td>" +
            "\n  <td style=\"text-align:center\"> &nbsp; --- </td>" +
            "\n  <td style=\"text-align:left\"> &nbsp; " + String(freq, format=
            "9.4f") + " </td>" + "\n  <td style=\"text-align:left\"> &nbsp; "
             + String(damp, format="9.4f") + " </td>\n</tr>", tempFile);

        end for;

        print("</table>\n", tempFile);
      end printTable;
    algorithm
      if size(c,1) < 1 then
         return;
      end if;

      // Sort complex vector
      if sort then
        (cSorted,cIndex) :=Modelica.ComplexMath.Vectors.sort(
              c,
              ascending,
              sortFrequency);
      else
        cSorted :=c;
        cIndex :=1:nc;
      end if;

      // Number of real zeros
      nReal := Modelica_LinearSystems2.Internal.numberOfRealZeros(cSorted);

      // Remove temporary file, if it exists
      Files.removeFile(tempFile);

      // Print heading
      // Following doesn't work in Dymola
      //Streams.print("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">", tempFile);
      print("<html>", tempFile);
      print("<style type=\"text/css\">", tempFile);
      print("* { font-size: 10pt; font-family: Arial,sans-serif; }", tempFile);
      print("</style>", tempFile);

      print("<table style=\"background-color:rgb(100, 100, 100); margin:20px 0 20px 20px;\" "
             + "cellpadding=\"3\" border=\"0\" cellspacing=\"1\">", tempFile);
      print("<caption>" + heading + "</caption>", tempFile);
      print("<tr style=\"background-color:rgb(230, 230, 230); text-align:center;\">" +
            "\n  <td> number </td>\n  <td>" + name + "</td>\n  <td> time constant [s] </td>" +
            "\n  <td> freq. [Hz] </td>\n  <td> damping </td>\n</tr>", tempFile);

      // Print values of complex vector
      printTable(cSorted, nReal, tempFile);

      // Print end of file
      print("</html>", tempFile);

      // Read file into output window and remove temporary file
      readFile(tempFile);
      Files.removeFile(tempFile);
    end printHTML;

    annotation (Documentation(info="<html>
<p>
This package provides functions operating on vectors of complex numbers.
</p>
</html>"));
  end Vectors;

  encapsulated package Matrices
    "Package of functions operating on matrices of complex numbers"
    extends Modelica.Icons.Package;
    import Modelica;
    import Modelica_LinearSystems2;

    function print "Print matrix"
      extends Modelica.Icons.Function;

      import Modelica;
      import Modelica.Utilities.Strings;
      import Complex;
      import Modelica_LinearSystems2.StateSpace;

      input Complex M[:,:];
      input Integer significantDigits=6
        "Number of significant digits that are shown";
      input String name="M" "Independent variable name used for printing";
      output String s="";
    protected
      String blanks=Strings.repeat(significantDigits);
      String space=Strings.repeat(8);
      String space2=Strings.repeat(3);
      Integer r=size(M, 1);
      Integer c=size(M, 2);

    algorithm
      if r == 0 or c == 0 then
        s := name + " = []";
      else
        s := "\n" + name + " = \n";
        for i in 1:r loop
          s := s + space;
          for j in 1:c loop
            if Modelica.ComplexMath.'abs'(M[i, j]) >= 0 then
              s := s + " ";
            end if;
            s :=s + String(M[i, j], significantDigits=significantDigits) + Strings.repeat(significantDigits + 8 - Strings.length(String(Modelica.ComplexMath.'abs'(M[i, j]))));

          end for;
          s := s + "\n";
        end for;

      end if;
    end print;

    encapsulated function matMatMul "Multiply two complex matrices"
      extends Modelica.Icons.Function;

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;

      input Complex m1[:,:] "Complex matrix 1";
      input Complex m2[size(m1, 2),:] "Complex matrix 2";
      output Complex m3[size(m1, 1),size(m2, 2)] "= m1*m2";
    //   Complex j=Complex.j();
    //   Real m1Re[size(m1,1),size(m1,2)]=Re(m1);
    //   Real m1Im[size(m1,1),size(m1,2)]=Im(m1);
    //   Real m2Re[size(m2,1),size(m2,2)]=Re(m2);
    //   Real m2Im[size(m2,1),size(m2,2)]=Im(m2);
    protected
      Integer l1;
      Integer l2;
    algorithm
    //  m3 := m1[:, :].re*m2[:, :].re - m1[:, :].im*m2[:, :].im + j*(m1[:, :].re*m2[:,:].im + m1[:, :].im*m2[:, :].re);
    //  m3 :=m1Re*m2Re - m1Im*m2Im + j*(m1Re*m2Im + m1Im*m2Re);

     for l1 in 1:size(m1,1) loop
       for l2 in 1:size(m2,2) loop
         m3[l1,l2] :=Complex.'*'.scalarProduct(m1[l1, :], m2[:, l2]);
       end for;
       end for;

    end matMatMul;

    encapsulated function matVecMul
      "Multiply a complex matrices with a complex vector"
      extends Modelica.Icons.Function;

      import Modelica;
      import Complex;
      import Modelica_LinearSystems2;

      input Complex m[:,:] "Complex matrix";
      input Complex vi[size(m, 2)] "Complex vector";
      output Complex vo[size(m, 1)] "= m*vi";
    //   Complex j=Complex.j();
    //   Real Rem[size(m, 1),size(m, 2)]=Re(m);
    //   Real Imm[size(m, 1),size(m, 2)]=Im(m);
    //   Real Revi[size(vi, 1)]=Re(vi);
    //   Real Imvi[size(vi, 1)]=Im(vi);
    protected
      Integer l1;

    algorithm
    //  vo := Rem*Revi - Imm*Imvi + j*(Rem*Imvi + Imm*Revi);

    //    for l1 in 1:size(m, 1) loop
    //      vo[l1] := Complex(0);
    //      for l2 in 1:size(m, 2) loop
    //        vo[l1] := vo[l1] + m[l1, l2]*vi[l2];
    //      end for;//l2
    //    end for;  //l1

    for l1 in 1:size(m, 1) loop
      vo[l1] :=Complex.'*'.scalarProduct(m[l1, :], vi);
    end for;

    end matVecMul;

    annotation (Documentation(info="<html>
<p>
This package provides functions operating on matrices of complex numbers.
</p>
</html>"));
  end Matrices;

  encapsulated function j "Returns sqrt(-1)"
    extends Modelica.Icons.Function;
    import Modelica;
    import Modelica_LinearSystems2.Math.Complex;

    output Complex c "= sqrt(-1)";
  algorithm
    c := Complex(0,1);
    annotation(Inline=true, Icon(graphics={
          Ellipse(
            lineColor={255,0,0},
            extent={{-100,-100},{100,100}},
            pattern=LinePattern.Dash,
            lineThickness=0.5)}));
  end j;

  encapsulated function eigenValues
    "Compute eingenvalues of a matrix A, using lapack routine dgeevx"
    extends Modelica.Icons.Function;

    import Modelica;
    import Complex;
    import Modelica_LinearSystems2.Math.Matrices.LAPACK;
    import Modelica_LinearSystems2.Math.Matrices.LAPACK.dgeevx;

    input Real A[:,size(A, 1)] "Real matrix";
    output Complex eigval[size(A, 1)]
      "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

  protected
    Integer nx=size(A, 1) "Number of states";

    Real alphaReal[nx];
    Real alphaImag[nx];

    Integer info;

  algorithm
   // Compute eigenvalues
    if size(A, 1) > 0 then
       (alphaReal,alphaImag,,,,info) := dgeevx(A);
       assert(info == 0,
         "Failed to compute eigenvalues with function eigenValues_dgeevx(..)");

       for i in 1:nx loop
         eigval[i].re := alphaReal[i];
         eigval[i].im := alphaImag[i];
       end for;
    end if;
    annotation (Documentation(info="<html>
<p>
Computes the invariant zeros of a system in state space form:
</p>
<pre>
   der(<b>x</b>) = <b>A</b>*<b>x</b> + <b>B</b>*<b>u</b>
        <b>y</b> = <b>C</b>*<b>x</b> + <b>D</b>*<b>u</b>
</pre>
<p>
The invariant zeros of this system are defined as the variables
z that make the following matrix singular:
</p>
<pre>
    | <b>A</b> <b>B</b> |     | <b>I</b> <b>0</b> |
    |     | - z*|     |
    | <b>C</b> <b>D</b> |     | <b>0</b> <b>0</b> |
</pre>
<p>
where <b>I</b> is the identity matrix of the same size as <b>A</b>
and <b>0</b> are zero matrices of appropriate dimensions.
</p>
<p>
Currently, there is the restriction that the number of
inputs and the number of outputs must be identical.
</p>
</html>"));
  end eigenValues;

  encapsulated function eigenVectors
    "Calculate the rigth eigenvectors of a linear state space system and write them columnwise in a matrix."
    extends Modelica.Icons.Function;

    import Modelica;
    import Modelica.Math.Matrices.LAPACK;
    import Complex;
    import Modelica.ComplexMath.j;

    input Real A[:,size(A, 1)] "real square matrix";
    output Complex eigvec[size(A, 1),size(A, 2)] "eigen values of the system";
    output Complex eigval[size(A, 1)]=fill(Complex(0), size(A, 1))
      "eigen values of the system";
  protected
    Integer info;
    Real eigvecRe[size(A, 1),size(A, 2)];
    Real eigvalRe[size(A, 1)]=fill(0, size(A, 1));
    Real eigvalIm[size(A, 1)]=fill(0, size(A, 1));
    Integer n=size(A, 1);
    Integer i;
  algorithm
    if size(A, 1) > 0 then

      (eigvalRe,eigvalIm,eigvecRe,info) := LAPACK.dgeev(A);
      for i in 1:size(A, 1) loop
        eigval[i].re := eigvalRe[i];
        eigval[i].im := eigvalIm[i];
      end for;

      assert(info == 0, "Calculating the eigen values with function
  \"StateSpace.Analysis.eigenVectors\" is not possible, since the
  numerical algorithm does not converge.");

      i := 1;
      while i <= n loop
        if abs(eigvalIm[i]) > 0 then
          for ii in 1:n loop
            eigvec[ii, i] := eigvecRe[ii, i] + j*eigvecRe[ii, i + 1];
            eigvec[ii, i + 1] := eigvecRe[ii, i] - j*eigvecRe[ii, i + 1];
          end for;
          i := i + 2;
        else
          for ii in 1:n loop
            eigvec[ii, i] := Complex(1)*eigvecRe[ii, i];
          end for;
          i := i + 1;
        end if;
      end while;

    end if;

    annotation (Documentation(info="<html>
  <h4>Syntax</h4>
  <table>
  <tr> <td align=right>  (eigenvectors, eigenvalues) </td><td align=center> =  </td>  <td> StateSpace.Analysis.<b>eigenVectors</b>(ss, onlyEigenvectors)  </td> </tr>
  </table>
  <h4>Description</h4>
  <p>
  Calculate the eigenvectors and optionally (onlyEigenvectors=false) the eigenvalues of a state space system. The output <tt>eigenvectors</tt> is a matrix with the same dimension as matrix <b>ss.A</b>. Just like in <a href=\"modelica://Modelica.Math.Matrices.eigenValues\">Modelica.Math.Matrices.eigenValues</a>, if the i-th eigenvalue has an imaginary part, then <tt>eigenvectors</tt>[:,i] is the real and <tt>eigenvectors</tt>[:,i+1] is the imaginary part of the eigenvector of the i-th eigenvalue.<br>
  The eigenvalues are returned as a complex vector <tt>eigenvalues</tt>.


  </p>

  <h4>Example</h4>
  <blockquote><pre>
     Modelica_LinearSystems2.StateSpace ss=Modelica_LinearSystems2.StateSpace(
        A=[-1,1;-1,-1],
        B=[1;1],
        C=[1,1],
        D=[0]);

     Real eigenvectors[2,2];
     Complex eigenvalues[2];

  <b>algorithm</b>
    (eigenvectors, eigenvalues) = Modelica_LinearSystems2.StateSpace.Analysis.eigenVectors(ss, true);
  // eigenvectors = [0.707, 0; 0, 0.707]
  // eigenvalues = {-1 + 1j, -1 - 1j}

            |0.707 |         | 0.707 |
  i.e. v1 = |      |,   v2 = |       |
            |0.707i|         |-0.707i|
  </pre></blockquote>


  </html>"));
  end eigenVectors;

  encapsulated function frequency
    "Frequency and damping of conjugated complex pole pair"
    extends Modelica.Icons.Function;

    import Modelica;
    import Complex;
    input Complex c "Complex number";
    output Modelica.SIunits.Frequency f "Frequency of c (= c.im in Hz)";
    output Real damping "Damping of c (= c.re/c.im)";

  protected
    Real abs_ev=(c.re^2 + c.im^2)^0.5;
  algorithm
    f := if abs(c.im) > 10*Modelica.Constants.eps then abs_ev/(2*Modelica.Constants.pi) else 0;
    damping := if abs(c.im) > 10*Modelica.Constants.eps then if abs_ev > Modelica.Constants.eps then -c.re/abs_ev else 0.0 else
      1.0;
    annotation(Inline=true);
  end frequency;

  encapsulated package Internal
    "Package of internal functions operating on complex number (for advanced users only)"
    extends Modelica.Icons.InternalPackage;
    import Modelica;

    function eigenValues_dhseqr
      "Compute eingenvalues of a upper Hessenberg matrix using lapack routine DHSEQR"
      extends Modelica.Icons.Function;

      import Modelica;
      import Complex;
      import
        Modelica_LinearSystems2.Math.Matrices.Internal.eigenvaluesHessenberg;

      input Real H[:,size(H, 1)] "Real upper Hessenberg matrix";
      output Complex zeros[size(H, 1)]
        "Finite, invariant zeros of ss; size(Zeros,1) <= size(ss.A,1)";

    protected
      Integer nx=size(H, 1) "Number of states";
      Real alphaReal[nx];
      Real alphaImag[nx];
      Integer info;

    algorithm
      (alphaReal,alphaImag,info) := eigenvaluesHessenberg(H);
      assert(info == 0,
        "Failed to compute eigenvalues with function Internal.eigenValues_dhseqr(..)");

      for i in 1:nx loop
        zeros[i].re := alphaReal[i];
        zeros[i].im := alphaImag[i];
      end for;

    end eigenValues_dhseqr;

    function C_transpose "Computes the transposed matrix of a complex matrix"
      extends Modelica.Icons.Function;

      import Complex;
      import Modelica.ComplexMath.j;

      input Complex C[:,:];
      output Complex CT[size(C, 2),size(C, 1)];
    protected
      Integer l1;
      Integer l2;
    algorithm
    //  CT := Complex(1)*transpose(C[:,:].re)-j*transpose(C[:,:].im);// too slow
      for l1 in 1:size(C, 1) loop
        for l2 in 1:size(C, 2) loop
    //      CT[l2, l1] := Re(C[l1, l2]) - j*Im(C[l1, l2]);
          CT[l2, l1] :=Modelica.ComplexMath.conj(C[l1, l2]);
        end for;
      end for;

    end C_transpose;

    function frobeniusNorm "Return the Frobenius norm of a matrix"
      extends Modelica.Icons.Function;

      import Complex;
      input Complex A[:,:] "Input matrix";
      output Real result=0.0 "frobenius norm of matrix A";
    algorithm
      for i1 in 1:size(A, 1) loop
        for i2 in 1:size(A, 2) loop
          result :=result + Modelica.ComplexMath.real(A[i1, i2]*Modelica.ComplexMath.conj(A[i1, i2]));
        end for;
      end for;
      result := sqrt(result);
    end frobeniusNorm;

    annotation (Documentation(info="<html>
<p>
Generally, the functions in this package should not be used by the user.
</p>
<p>
This package contains functions which cannot be used in an arbitrary
way and require particular knowledge.
Therefore, only advanced users should deal with such a functions.
</p>
</html>"));
  end Internal;

  annotation (
    Documentation(info="<html>
<p>
This package contains some <b>utility functions</b>
operating on complex numbers (such as frequency(..)), as well as
functions operating on vectors and matrices of complex numbers.
For general information about the complex numbers usage, please refer to
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.GettingStarted.ComplexNumbers\">GettingStarted.ComplexNumbers</a>.
</p>
<p>
Example (note: \"j\" in the comments is defined as j=sqrt(-1)):
</p>
<blockquote>
<pre>
  <span style=\"font-family:'Courier New,courier'; color:#0000ff;\">import </span><span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Modelica.Utilities.Streams.print</span>;
  <span style=\"font-family:'Courier New,courier'; color:#0000ff;\">import </span>LS2 = <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Modelica_LinearSystems2</span>;
  <span style=\"font-family:'Courier New,courier'; color:#0000ff;\">import </span><span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span>;

  <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span> c1 = <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span>(re=2, im=3) <span style=\"font-family:'Courier New,courier'; color:#006400;\">\"= 2 + 3j\"</span>;
  <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span> c2 = <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span>(3,4) <span style=\"font-family:'Courier New,courier'; color:#006400;\">\"= 3 + 4j\"</span>;
  <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span> c3;

<span style=\"font-family:'Courier New,courier'; color:#0000ff;\">algorithm </span>
  c3 := c1 + c2;
  <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">print</span>(<span style=\"font-family:'Courier New,courier'; color:#006400;\">\"c3 = \"</span> + <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">LS2.Math.Complex.'String'</span>(c3));
  <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">print</span>(<span style=\"font-family:'Courier New,courier'; color:#006400;\">\"c3 = \"</span> + <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">LS2.Math.Complex.'String'</span>(c3,<span style=\"font-family:'Courier New,courier'; color:#006400;\">\"i\"</span>));

  <span style=\"color:#006400;\">// This gives the following print-out:</span>
  <span style=\"color:#006400;\">// c3 = 5 + 7j</span>
  <span style=\"color:#006400;\">// c3 = 5 + 7i</span>
  ...
<span style=\"font-family:'Courier New,courier'; color:#0000ff;\">end </span>functionName;
</pre>
</blockquote>

<p>
The utility functions are written in such a way that it
is convenient to work with them, once operator overloading
is provided in Modelica and Modelica tools. Example:
</p>
<blockquote>
<pre>
  <span style=\"color:#006400;\">// Assume operator overloading is available</span>:
  <span style=\"font-family:'Courier New,courier'; color:#0000ff;\">import</span> <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span>;
  <span style=\"font-family:'Courier New,courier'; color:#0000ff;\">import</span> <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Modelica.ComplexMath.j</span>;

  <span style=\"font-family:'Courier New,courier'; color:#ff0000;\">Complex</span> c4 = -2 + 5*j;
</pre>
</blockquote>
<p>
The utility functions are implemented so that they can be easily
inlined by a tool. In such a case, the above statement will
not lead to any overhead.
</p>
</html>"));
end Complex;
