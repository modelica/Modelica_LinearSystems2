within Modelica_LinearSystems2;
record DiscreteZerosAndPoles
  "Discrete zeros and poles description of a single input, single output system (data + operations)"
  extends Modelica.Icons.Record;

  Real k=1.0 "Multiplicative factor of transfer function" 
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real n1[:] "[p^0] coefficients of 1st order numerator polynomials" 
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real n2[:,2] "[p,p^0] coefficients of 2nd order numerator polynomials" 
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real d1[:] "[p^0] coefficients of 1st order denominator polynomials" 
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
  Real d2[:,2] "[p,p^0] coefficients of 2nd order denominator polynomials" 
      annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));

    Modelica.SIunits.Time Ts "Sample time" 
       annotation(Dialog(group="Data used to construct discrete from continuous system"));

  Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
    "Discretization method" 
        annotation(Dialog(group="Data used to construct discrete from continuous system"));

/* If the numerator polynomial has no coefficients, the transfer function
   is zero. The denominator polynomial must always have at
   least one coefficient, such as {1}
*/

  String uName="u" "Name of input signal"    annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal"  annotation(Dialog(group="Signal names"));

  encapsulated operator 'constructor' "Generate a ZerosAndPoles object"
    import Modelica_LinearSystems2;

    encapsulated function fromReal
      "Generate a DiscreteZerosAndPoles transfer function from a Real value"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

      input Real r "Value of Real variable";
      input Modelica.SIunits.Time Ts=1 "Sample time";
         input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteZerosAndPoles dzp(
        redeclare Real n1[0],
        redeclare Real n2[0,2],
        redeclare Real d1[0],
        redeclare Real d2[0,2]) "= r";

    algorithm
      dzp.k := r;
      dzp.Ts := Ts;
      dzp.method := method;
      dzp.uName := uName;
      dzp.yName := yName;

      annotation (overloadsConstructor=true);
    end fromReal;

  function fromZerosAndPoles
      "Generate a DiscreteZerosAndPoles data record from a continuous ZerosAndPoles represenation"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteStateSpace;

    input ZerosAndPoles zp "continuous zeros and poles transfer function";
    input Modelica.SIunits.Time Ts "Sample time" 
         annotation(Dialog(group="Data used to construct discrete from continuous system"));

    input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method" 
          annotation(Dialog(group="Data used to construct discrete from continuous system"));

    output DiscreteZerosAndPoles dzp;
    protected
    StateSpace ss = StateSpace(zp);
    Modelica_LinearSystems2.DiscreteStateSpace dss=
                             Modelica_LinearSystems2.DiscreteStateSpace(ss,Ts,method);

  algorithm
    dzp := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss);
    annotation (overloadsConstructor=true);
  end fromZerosAndPoles;

  encapsulated function fromPolesAndZeros
      "Generate a DiscreteZerosAndPoles transfer function from a set of zeros and poles"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica.Utilities.Streams.print;

    input Complex z[:]=fill(Modelica_LinearSystems2.Math.Complex(0), 0)
        "Zeros (Complex vector of numerator zeros)";
    input Complex p[:]=fill(Modelica_LinearSystems2.Math.Complex(0), 0)
        "Poles (Complex vector of denominator zeros)";
    input Real k=1.0 "Constant multiplied with transfer function";
    input Modelica.SIunits.Time Ts "Sample time";
       input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
    input String uName="" "input name";
    input String yName="" "output name";
    output DiscreteZerosAndPoles dzp(
      redeclare Real n1[Internal.numberOfRealZeros(z)],
      redeclare Real n2[integer((size(z, 1) - Internal.numberOfRealZeros(z))/2),
        2],
      redeclare Real d1[Internal.numberOfRealZeros(p)],
      redeclare Real d2[integer((size(p, 1) - Internal.numberOfRealZeros(p))/2),
        2]) "ZerosAndPoles transfer functions of the zeros, poles and k";

    protected
    Integer n_n1=size(dzp.n1, 1);
    Integer n_d1=size(dzp.d1, 1);
    Integer n_n2=size(dzp.n2, 1);
    Integer n_d2=size(dzp.d2, 1);
    Integer i;
    Integer j;
    Complex z_reordered[size(z, 1)] "Reordered zeros";
    Complex p_reordered[size(p, 1)] "Reordered poles";
    Integer nz_real "Number of real zeros";
    Integer np_real "Number of real poles";
  algorithm
    //Extra input (string) added
    (z_reordered,nz_real) := Internal.reorderZeros(z, "");
    (p_reordered,np_real) := Internal.reorderZeros(p, "");

    // Numerator
    for i in 1:n_n1 loop
      dzp.n1[i] := -z_reordered[i].re;
    end for;

    j := 1;
    for i in n_n1 + 1:2:size(z, 1) loop
      dzp.n2[j, :] := {-2*z_reordered[i].re,z_reordered[i].re^2 + z_reordered[i].im
        ^2};
      j := j + 1;
    end for;

    // Denominator
    for i in 1:n_d1 loop
      dzp.d1[i] := -p_reordered[i].re;
    end for;

    j := 1;
    for i in n_d1 + 1:2:size(p, 1) loop
      dzp.d2[j, :] := {-2*p_reordered[i].re,p_reordered[i].re^2 + p_reordered[i].im
        ^2};
      j := j + 1;
    end for;

    dzp.Ts := Ts;
    dzp.method := method;
    dzp.k := k;
    dzp.uName := uName;
    dzp.yName := yName;

    annotation (Documentation(info="<html>

</html>"));
  end fromPolesAndZeros;

    function fromDiscreteTransferFunction = 
        Modelica_LinearSystems2.DiscreteTransferFunction.Conversion.toDiscreteZerosAndPoles;
    encapsulated function fromFactorization
      "Generate a DiscreteZerosAndPoles object from first and second order polynomials"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

      input Real n1[:]=fill(0, 0)
        "[p^0] coefficients of 1st order numerator polynomials" 
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real n2[:,2]=fill(
              0,
              0,
              2) "[p,p^0] coefficients of 2nd order numerator polynomials" 
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real d1[:]=fill(0, 0)
        "[p^0] coefficients of 1st order denominator polynomials" 
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real d2[:,2]=fill(
              0,
              0,
              2) "[p,p^0] coefficients of 2nd order denominator polynomials" 
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Real k=1.0 "Multiplicative factor of transfer function" 
           annotation(Dialog(group="y = k*(product(p+n1[i]) * product(p^2+n2[i,1]*p+n2[i,2])) / (product(p+d1[i])*product(p^2+d2[i,1]*p+d2[i,2])) *u"));
      input Modelica.SIunits.Time Ts=1 "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteZerosAndPoles dzp(
        redeclare Real n1[size(n1, 1)],
        redeclare Real n2[size(n2, 1),2],
        redeclare Real d1[size(d1, 1)],
        redeclare Real d2[size(d2, 1),2]) "ZerosAndPoles transfer function";
    algorithm
      dzp.n1 := n1;
      dzp.n2 := n2;
      dzp.d1 := d1;
      dzp.d2 := d2;
      dzp.k := k;
      dzp.Ts := Ts;
      dzp.method := method;
      dzp.uName := uName;
      dzp.yName := yName;
    end fromFactorization;

  end 'constructor';

encapsulated operator '-'
  function subtract "Subtract two DiscreteZerosAndPoles (dzp1 - dzp2)"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.Math.Complex;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;

    protected
    Integer size_z1n1=size(dzp1.n1, 1);
    Integer size_z1d1=size(dzp1.d1, 1);
    Integer size_z1n2=size(dzp1.n2, 1);
    Integer size_z1d2=size(dzp1.d2, 1);
    Integer size_z2n1=size(dzp2.n1, 1);
    Integer size_z2d1=size(dzp2.d1, 1);
    Integer size_z2n2=size(dzp2.n2, 1);
    Integer size_z2d2=size(dzp2.d2, 1);
    Polynomial p1;
    Polynomial p2;
    Polynomial p3;
    Complex numZeros[:];
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 +
        size_z2d2));
    Real k;

    output DiscreteZerosAndPoles result "= dzp1/dzp2";
  algorithm
    assert(abs(dzp1.Ts-dzp2.Ts)<=Modelica.Constants.eps,"Two discrete zeros-and-poles systems must have the same sample time Ts for subtraction with \"-.subtract\".");
    result.Ts := dzp1.Ts;
     if dzp1 == dzp2 then
       result := DiscreteZerosAndPoles(0);
     else
      p1 := Polynomial(1);
      p2 := Polynomial(1);

      for i in 1:size_z1n1 loop
        p1 := p1*Polynomial({1,dzp1.n1[i]});
      end for;
      for i in 1:size_z1n2 loop
        p1 := p1*Polynomial(cat(1, {1}, dzp1.n2[i, :]));
      end for;
      for i in 1:size_z2d1 loop
        p1 := p1*Polynomial({1,dzp2.d1[i]});
      end for;
      for i in 1:size_z2d2 loop
        p1 := p1*Polynomial(cat(1, {1}, dzp2.d2[i, :]));
      end for;

      for i in 1:size_z2n1 loop
        p2 := p2*Polynomial({1,dzp2.n1[i]});
      end for;
      for i in 1:size_z2n2 loop
        p2 := p2*Polynomial(cat(1, {1}, dzp2.n2[i, :]));
      end for;
      for i in 1:size_z1d1 loop
        p2 := p2*Polynomial({1,dzp1.d1[i]});
      end for;
      for i in 1:size_z1d2 loop
        p2 := p2*Polynomial(cat(1, {1}, dzp1.d2[i, :]));
      end for;

      p3 := dzp1.k*p1 - dzp2.k*p2;
      k := 0;
      for i in size(p3.c, 1):-1:1 loop
        if abs(p3.c[i]) > Modelica.Constants.eps then
          k := p3.c[i];
        end if;
      end for;
      numZeros := Polynomial.roots(p3);
      result := DiscreteZerosAndPoles(numZeros, dummy, k, dzp1.Ts, dzp1.method);

      result.d1 := cat(1, dzp1.d1, dzp2.d1);
      result.d2 := cat(1, dzp1.d2, dzp2.d2);

    end if;
  end subtract;

  function negate "Unary minus (multiply transfer function by -1)"
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

     input DiscreteZerosAndPoles dzp;
     output DiscreteZerosAndPoles result(n1=dzp.n1, n2=dzp.n2, d1=dzp.d1, d2=dzp.d2, k=-dzp.k, Ts=dzp.Ts, method=dzp.method) "= -dzp";
  algorithm
  end negate;
end '-';

encapsulated operator '*'
function 'dzp*dzp'
      "Multiply two DiscreteZerosAndPoles transfer functions (dzp1 * dzp2)"

      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2;

  input DiscreteZerosAndPoles dzp1;
  input DiscreteZerosAndPoles dzp2;
  input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";
  input String uName=dzp1.uName "input name";
  input String yName=dzp2.yName "output name";

  output DiscreteZerosAndPoles result "= dzp1 * dzp2";
algorithm
  assert(abs(dzp1.Ts-dzp2.Ts)<=Modelica.Constants.eps,"Two discrete zeros-and-poles systems must have the same sample time Ts for multiplication with \"*\".");
  result.Ts := dzp1.Ts;
  if dzp1 == DiscreteZerosAndPoles(0) or dzp2 == DiscreteZerosAndPoles(0) then
    result := DiscreteZerosAndPoles(0);
  else
    result.n1 := cat(
      1,
      dzp1.n1,
      dzp2.n1);
    result.n2 := cat(
      1,
      dzp1.n2,
      dzp2.n2);
    result.d1 := cat(
      1,
      dzp1.d1,
      dzp2.d1);
    result.d2 := cat(
      1,
      dzp1.d2,
      dzp2.d2);
    result.k := dzp1.k*dzp2.k;
  end if;
    result.Ts := dzp1.Ts;
    result.method := dzp1.method;
    result.uName := uName;
    result.yName := yName;

end 'dzp*dzp';

function 'r*dzp'
      "Multiply a real number with a discrete DiscreteZerosAndPoles transfer function  (r * dzp2)"

      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2;

  input Real r;
  input DiscreteZerosAndPoles dzp;
  input String uName=dzp.uName "input name";
  input String yName=dzp.yName "output name";

  output DiscreteZerosAndPoles result "= r * dzp1";
algorithm
  if r == 0 or dzp == DiscreteZerosAndPoles(0) then
    result := DiscreteZerosAndPoles(0);
  else
    result := dzp;
    result.k := r*dzp.k;
  end if;
    result.Ts := dzp.Ts;
    result.method := dzp.method;
    result.uName := uName;
    result.yName := yName;

end 'r*dzp';
end '*';


  encapsulated operator function '/'
    "Divide two discrete transfer functions in zeros and poles representation (dzp1 / dzp2)"
    import Modelica;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;

    output DiscreteZerosAndPoles result "= dzp1/dzp2";

  algorithm
    assert(abs(dzp2.k) > 100*Modelica.Constants.small, "dzp2 in operator \"Modelica_LinearSystems2.TransferFunction.'/'()\" may not be zero");
    if max({size(dzp1.n1, 1),size(dzp1.n2, 1),size(dzp1.d1, 1),size(dzp1.d2, 1)}) >0 and 
      max({size(dzp2.n1, 1),size(dzp2.n2, 1),size(dzp2.d1, 1),size(dzp2.d2, 1)}) > 0 then
        assert(abs(dzp1.Ts - dzp2.Ts) <= Modelica.Constants.eps, "Two discrete zeros-and-poles systems must have the same sample time Ts for division with \"/\".");
    end if;
     result.Ts := dzp1.Ts;
    result.method := dzp1.method;

    if dzp1 == DiscreteZerosAndPoles(0) then
      result := DiscreteZerosAndPoles(0);
    else
      result.n1 := cat(
        1,
        dzp1.n1,
        dzp2.d1);
      result.n2 := cat(
        1,
        dzp1.n2,
        dzp2.d2);
      result.d1 := cat(
        1,
        dzp1.d1,
        dzp2.n1);
      result.d2 := cat(
        1,
        dzp1.d2,
        dzp2.n2);
      result.k := dzp1.k/dzp2.k;
    end if;
  end '/';

  encapsulated operator function '+'
    "Addition of to discrete transfer functions dzp1 + dzp2, i.e. parallel connection of two transfer functions (= inputs are the same, outputs of the two systems are added)"

    import Modelica;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.Math.Complex;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;

  protected
    Integer size_z1n1=size(dzp1.n1, 1);
    Integer size_z1d1=size(dzp1.d1, 1);
    Integer size_z1n2=size(dzp1.n2, 1);
    Integer size_z1d2=size(dzp1.d2, 1);
    Integer size_z2n1=size(dzp2.n1, 1);
    Integer size_z2d1=size(dzp2.d1, 1);
    Integer size_z2n2=size(dzp2.n2, 1);
    Integer size_z2d2=size(dzp2.d2, 1);
    Polynomial p1;
    Polynomial p2;
    Polynomial p3;
    Complex numZeros[:];
    Complex dummy[:]=fill(Complex(1), size_z1d1 + size_z2d1 + 2*(size_z1d2 +
        size_z2d2));
    Real k;

    output DiscreteZerosAndPoles result "= dzp1+dzp2";

  algorithm
    if max({size(dzp1.n1,1),size(dzp1.n2,1),size(dzp1.d1,1),size(dzp1.d2,1)})>0 and max({size(dzp2.n1,1),size(dzp2.n2,1),size(dzp2.d1,1),size(dzp2.d2,1)})>0 then
      assert(abs(dzp1.Ts-dzp2.Ts)<=Modelica.Constants.eps,"Two discrete zeros-and-poles systems must have the same sample time Ts for addition with \"+\".");
    end if;
    result.Ts := dzp1.Ts;

    if dzp1 == -dzp2 then
      result := DiscreteZerosAndPoles(0,Ts=  dzp1.Ts);
    else

      p1 := Polynomial(1);
      p2 := Polynomial(1);

      for i in 1:size_z1n1 loop
        p1 := p1*Polynomial({1,dzp1.n1[i]});
      end for;
      for i in 1:size_z1n2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          dzp1.n2[i, :]));
      end for;
      for i in 1:size_z2d1 loop
        p1 := p1*Polynomial({1,dzp2.d1[i]});
      end for;
      for i in 1:size_z2d2 loop
        p1 := p1*Polynomial(cat(
          1,
          {1},
          dzp2.d2[i, :]));
      end for;

      for i in 1:size_z2n1 loop
        p2 := p2*Polynomial({1,dzp2.n1[i]});
      end for;
      for i in 1:size_z2n2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          dzp2.n2[i, :]));
      end for;
      for i in 1:size_z1d1 loop
        p2 := p2*Polynomial({1,dzp1.d1[i]});
      end for;
      for i in 1:size_z1d2 loop
        p2 := p2*Polynomial(cat(
          1,
          {1},
          dzp1.d2[i, :]));
      end for;

      p3 := dzp1.k*p1 + dzp2.k*p2;
      k := p3.c[1];
      numZeros := Polynomial.roots(p3);
      result := DiscreteZerosAndPoles(numZeros, dummy, k, dzp1.Ts, dzp1.method);

      result.d1 := cat(
        1,
        dzp1.d1,
        dzp2.d1);
      result.d2 := cat(
        1,
        dzp1.d2,
        dzp2.d2);
    end if;

  end '+';

  encapsulated operator function '^'
    "Integer power of DiscreteZerosAndPoles (dzp^k)"

    import Modelica;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp;
    input Integer k;

    output DiscreteZerosAndPoles result(
      redeclare Real n1[k*size(dzp.n1, 1)],
      redeclare Real n2[k*size(dzp.n2, 1),2],
      redeclare Real d1[k*size(dzp.d1, 1)],
      redeclare Real d2[k*size(dzp.d2, 1),2]) "= dzp^k";
  protected
    Integer size_n1=size(dzp.n1, 1);
    Integer size_d1=size(dzp.d1, 1);
    Integer size_n2=size(dzp.n2, 1);
    Integer size_d2=size(dzp.d2, 1);
  algorithm
    result.Ts := dzp.Ts;
    result.method := dzp.method;
    for i in 1:k loop
      result.n1[(i - 1)*size_n1 + 1:i*size_n1] := dzp.n1;
      result.d1[(i - 1)*size_d1 + 1:i*size_d1] := dzp.d1;
      result.n2[(i - 1)*size_n2 + 1:i*size_n2, :] := dzp.n2;
      result.d2[(i - 1)*size_d2 + 1:i*size_d2, :] := dzp.d2;
    end for;
    result.k := dzp.k^k;

  end '^';

encapsulated operator function '=='
    "Check whether two DiscreteZerosAndPoles transfer functions are identical"
    import Modelica;
    import Modelica.Math;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp1;
    input DiscreteZerosAndPoles dzp2;
    input Real eps(min=0) = 0
      "Two numbers n1 and n2 are identical if abs(n1-n2) <= eps";

    output Boolean result "= dzp1 == dzp2";
algorithm
    result := Math.Vectors.isEqual(dzp1.n1,dzp2.n1,eps) and Math.Vectors.isEqual(dzp1.d1,dzp2.d1,eps) and Math.Matrices.isEqual(dzp1.n2,dzp2.n2,eps) and Math.Matrices.isEqual(dzp1.d2,dzp2.d2,eps) and (dzp1.k==dzp2.k) and abs(dzp1.Ts-dzp2.Ts)<=eps;

end '==';

  encapsulated operator function 'String'
    "Transform DiscreteZerosAndPoles transfer function into a String representation"
    import Modelica;
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;
    import Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Types.Method;

      input DiscreteZerosAndPoles dzp
      "DiscreteZerosAndPoles transfer function to be transformed in a String representation";
      input Integer significantDigits=6
      "Number of significant digits that are shown";
      input String name="q" "Independent variable name used for printing";
      output String s="";
  protected
      Integer num_order=size(dzp.n1, 1) + 2*size(dzp.n2, 1);
      Integer den_order=size(dzp.d1, 1) + 2*size(dzp.d2, 1);
  algorithm
      if num_order == 0 and den_order == 0 then
        s := String(dzp.k);

      else
         // construct string for multiplicative factor
        if dzp.k <> 1.0 or dzp.k == 1.0 and num_order == 0 then
          s := String(dzp.k);
          if num_order <> 0 then
            s := s + "*";
          end if;
        end if;

        if num_order <> 0 then
            // construct numerator string
          s := s + ZerosAndPoles.Internal.firstOrderToString(
                dzp.n1,
                significantDigits,
                name);
          if size(dzp.n2, 1) <> 0 then
            s := if size(dzp.n1, 1) > 0 then s + "*" +
              ZerosAndPoles.Internal.secondOrderToString(
                  dzp.n2,
                  significantDigits,
                  name) else s + ZerosAndPoles.Internal.secondOrderToString(
                  dzp.n2,
                  significantDigits,
                  name);
          end if;
        end if;

        if den_order <> 0 then
            // construct denominator string
          s := s + "/";
          if den_order > 1 then
            s := s + "( ";
          end if;
          s := s + ZerosAndPoles.Internal.firstOrderToString(
                dzp.d1,
                significantDigits,
                name);
          if size(dzp.d2, 1) <> 0 then
            if size(dzp.d1, 1) > 0 then
              s := s + "*";
            end if;
            s := s + ZerosAndPoles.Internal.secondOrderToString(
                  dzp.d2,
                  significantDigits,
                  name);
          end if;
          if den_order > 1 then
            s := s + " )";
          end if;
        end if;
          end if;
          s := s +"\n    Ts = " + String(dzp.Ts) + "\n    method ="+ Modelica_LinearSystems2.Internal.methodString(dzp.method);
    //    end toString;
  end 'String';

  encapsulated function q
    "Generate the DiscreteZerosAndPoles transfer function q"
    import Modelica;
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    output DiscreteZerosAndPoles dzp(
      redeclare Real n1[1],
      redeclare Real n2[0,2],
      redeclare Real d1[0],
      redeclare Real d2[0,2]);
  algorithm
    dzp.n1[1] := 0;
    annotation (Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  q </td><td align=center> =  </td>  <td>DiscreteZerosAndPoles.<b>q</b>()  </td> </tr>
 
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Generate the complex Laplace variable q=exp(s*T) as a DiscreteZerosAndPoles transfer function. It can be used for generating like 
<blockquote><pre>
        DiscreteZerosAndPoles dzp = q/(q^2 + q + 1)/(q + 1)
</pre></blockquote>

</p>
 

 
 
</html> "));
  end q;

  encapsulated package Analysis

  encapsulated function timeResponse
      "Calculate the time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Types.TimeResponse;

      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;  // Input/Output declarations of discrete time response functions
      input Modelica_LinearSystems2.Types.TimeResponse response=Modelica_LinearSystems2.Types.TimeResponse.Step;
      input Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";

    protected
      DiscreteStateSpace dss=DiscreteStateSpace(dzp);
      Real tSpanVar;

  algorithm
    // set sample time
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(dss);
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) := DiscreteStateSpace.Analysis.timeResponse(
          dss=dss,
          tSpan=tSpanVar,
          response=response,
          x0=x0);

      annotation (Documentation(info="<html>

</html>"));
  end timeResponse;

  encapsulated function impulseResponse
      "Calculate the impulse time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
      Real tSpanVar;
  algorithm

  // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dzp));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) := DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Types.TimeResponse.Impulse,
          x0=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)));

      annotation (interactive=true, Documentation(info="<html>

 
</html> "));
  end impulseResponse;

  encapsulated function stepResponse
      "Calculate the step time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
    extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
    Real tSpanVar;
  algorithm

  // set simulation time span
    if tSpan == 0 then
      tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(DiscreteStateSpace(dzp));
    else
      tSpanVar := tSpan;
    end if;

    (y,t,x_discrete) := DiscreteZerosAndPoles.Analysis.timeResponse(
      dzp=dzp,
      tSpan=tSpanVar,
      response=Modelica_LinearSystems2.Types.TimeResponse.Step,
      x0=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)));

    annotation (interactive=true, Documentation(info="<html>

</html> "));
  end stepResponse;

  encapsulated function rampResponse
      "Calculate the ramp time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
      Real tSpanVar;
  algorithm

  // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dzp));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) := DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Types.TimeResponse.Ramp,
          x0=zeros(DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)));

      annotation (interactive=true, Documentation(info="<html>

</html> "));
  end rampResponse;

  encapsulated function initialResponse
      "Calculate the initial time response of a discrete zeros-and-poles transfer function"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input Real x0[:]=fill(0, 0) "Initial state vector";

      // Input/Output declarations of time response functions:
      extends Modelica_LinearSystems2.Internal.timeResponseMask_zp_discrete;

    protected
      Real tSpanVar;
  algorithm

  // set simulation time span
      if tSpan == 0 then
        tSpanVar := DiscreteStateSpace.Internal.timeResponseSamples(
          DiscreteStateSpace(dzp));
      else
        tSpanVar := tSpan;
      end if;

      (y,t,x_discrete) := DiscreteZerosAndPoles.Analysis.timeResponse(
          dzp=dzp,
          tSpan=tSpanVar,
          response=Modelica_LinearSystems2.Types.TimeResponse.Initial,
          x0=x0);

      annotation (interactive=true, Documentation(info="<html>

 
 
</html> "));
  end initialResponse;

  encapsulated function denominatorDegree
      "Return denominator degree of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output Integer result;
  algorithm
    result := size(dzp.d1, 1) + 2*size(dzp.d2, 1);
    annotation (Documentation(info="<html>


</html> "));
  end denominatorDegree;

  encapsulated function numeratorDegree
      "Return numerator degree of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output Integer result;
  algorithm
    result := size(dzp.n1, 1) + 2*size(dzp.n2, 1);
    annotation (Documentation(info="<html>


</html> "));
  end numeratorDegree;

  encapsulated function evaluate
      "Evaluate a DiscreteZerosAndPoles transfer function at a given value of q"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.ZerosAndPoles.Internal;
      import Modelica_LinearSystems2.Math.Complex;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    input Complex q=Complex(0) "Complex value q where dzp is evaluated";
    input Real den_min=0 "|denominator(p)| is limited by den_min";
    output Complex y "= zp(p)";
    protected
    Complex j = Modelica_LinearSystems2.Math.Complex.j();
    Complex num;
    Complex den;
    Real abs_den;
  algorithm
    // Build numerator
    num := dzp.k+0*j;
    for i in 1:size(dzp.n1, 1) loop
       num := num*Internal.'p+a'(q, dzp.n1[i]);
    end for;
    for i in 1:size(dzp.n2, 1) loop
       num := num*Internal.'p^2+k[1]*p+k[2]'(q, dzp.n2[i, :]);
    end for;

    // Build denominator
    den := 1+0*j;
    for i in 1:size(dzp.d1, 1) loop
      den := den*Internal.'p+a'(q, dzp.d1[i]);
    end for;
    for i in 1:size(dzp.d2, 1) loop
      den := den*Internal.'p^2+k[1]*p+k[2]'(q, dzp.d2[i, :]);
    end for;

    // Build value of transfer function
    abs_den := Complex.'abs'(den);
    den := if abs_den >= den_min then den else -abs_den+0*j;
    y := num/den;
    annotation (Documentation(info="<html>

</html> "));
  end evaluate;

  end Analysis;

  encapsulated package Design
  end Design;

  encapsulated package Plot

  encapsulated function bode
      "Plot discrete zeros a-and-poles transfer function as bode plot"
      import Modelica;
      import Modelica.Utilities.Strings;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Math.Complex;

      import Modelica_LinearSystems2.Utilities.Plot;

      import SI = Modelica.SIunits;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles function to be plotted";
    input Integer nPoints(min=2) = 200 "Number of points";
    input Boolean autoRange=true
        "= true, if abszissa range is automatically determined";
    input Modelica.SIunits.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false" 
                                                      annotation(Dialog(enable=not autoRange));
    input Modelica.SIunits.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false"                                                annotation(Dialog(enable=not autoRange));

    input Boolean magnitude=true "= true, to plot the magnitude of dzp" 
                                                                       annotation(choices(__Dymola_checkBox=true));
    input Boolean phase=true "= true, to plot the pase of dzp" annotation(choices(__Dymola_checkBox=true));

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot of  dzp = "
           + String(dzp)));

    protected
    SI.AngularVelocity w[nPoints];
    Complex z[nPoints];
    SI.Frequency f[nPoints];
    Modelica.SIunits.Conversions.NonSIunits.Angle_deg phi[nPoints];
    Real A[nPoints];
    Boolean OK;
    Complex c;
    Modelica.SIunits.Angle phi_old;
    Complex numZeros[:];
    Complex denZeros[:];
    Complex numZerosZ[:];
    Complex denZerosZ[:];
    ZerosAndPoles zp=ZerosAndPoles(k=dzp.k, n1=dzp.n1, n2=dzp.n2, d1=dzp.d1, d2=dzp.d2);

    Plot.Records.Curve curves[2];
    Integer i;
    Plot.Records.Diagram diagram2[2];

  algorithm
          // Determine frequency vector f
    if autoRange then
       (numZerosZ,denZerosZ) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
    else
      numZerosZ := fill(Complex(0), 0);
      denZerosZ := fill(Complex(0), 0);
    end if;

    numZeros := fill(Complex(0),0);
    denZeros := fill(Complex(0),size(denZerosZ,1));
    for i in 1:size(denZerosZ,1) loop
      denZeros[i] := Complex.log(denZerosZ[i])/dzp.Ts;
    end for;

    f := Internal.frequencyVector(
          nPoints,
          autoRange,
          f_min,
          f_max,
          numZeros,
          denZeros);

    // Compute magnitude/phase at the frequency points
    phi_old := 0.0;
    for i in 1:nPoints loop
      w[i] := Modelica.SIunits.Conversions.from_Hz(f[i]);
      z[i] := Complex.exp(Complex(0,w[i]*dzp.Ts));
      c := ZerosAndPoles.Analysis.evaluate(
            zp,
            z[i],
            1e-10);
      A[i] := Complex.'abs'(c);
      phi_old := Complex.arg(c, phi_old);
      phi[i] := Modelica.SIunits.Conversions.to_deg(phi_old);

    end for;

   // Plot computed frequency response
   diagram2 := fill(defaultDiagram, 2);
    i := 0;
    if magnitude then
      i := i + 1;
      curves[i] := Plot.Records.Curve(
            x=f,
            y=A,
            autoLine=true);
      diagram2[i].curve := {curves[i]};
      diagram2[i].yLabel := "magnitude";
      if phase then
         diagram2[i].xLabel:="";
      end if;
    end if;

    if phase then
      i := i + 1;
      curves[i] := Plot.Records.Curve(
            x=f,
            y=phi,
            autoLine=true);
      diagram2[i].curve := {curves[i]};
      diagram2[i].yLabel := "phase [deg]";
      diagram2[i].logY := false;
      if magnitude then
        diagram2[i].heading:="";
     end if;
    end if;

      if magnitude and phase then
        Plot.diagramVector(diagram2, device);
      else
        Plot.diagram(diagram2[1], device);
      end if;

    annotation (interactive=true, Documentation(info="<html>



</html> "));
  end bode;

  encapsulated function timeResponse
      "Plot the time response of a system represented by a discrete zeros-and-poles transfer function. The response type is selectable"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.DiscreteZerosAndPoles dzp;
  //  input Real dt=0 "Sample time [s]";
    input Real tSpan=0 "Simulation time span [s]";

    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Step "Type of time response";
    input Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Time response of  dzp = "
           + String(dzp)));

    protected
    Plot.Records.Curve curve;
    Plot.Records.Diagram diagram2;
    Real y[:,1,1] "Output response";
    Real t[:] "Time vector: (number of samples)";

    Real yy[:] "Output response";
    Real tt[:] "Time vector: (number of samples)";

  algorithm
    (y,t) := DiscreteZerosAndPoles.Analysis.timeResponse(
      dzp,
      tSpan,
      response,
      x0);

    tt := fill(0,2*size(t,1)-1);
    yy := fill(0,2*size(t,1)-1);

    for i in 1:size(t,1)-1 loop
      tt[2*i-1] := t[i];
      tt[2*i] := t[i+1];
      yy[2*i-1] := y[i,1,1];
      yy[2*i] := y[i,1,1];
    end for;
    tt[size(tt,1)] := t[size(t,1)];
    yy[size(tt,1)] := y[size(t,1),1,1];

    curve := Plot.Records.Curve(
      x=tt,
      y=yy,
      legend="y",
      autoLine=true);
    diagram2 := defaultDiagram;
    diagram2.curve := {curve};

    Plot.diagram(diagram2, device);
    annotation (interactive=true, Documentation(info="<html>


</html> "));
  end timeResponse;

  encapsulated function impulse
      "Impulse response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

      import Modelica_LinearSystems2.Utilities.Plot;

      input DiscreteZerosAndPoles dzp "zeros-and-poles transfer function";
      input Real tSpan=0 "Simulation time span [s]";

      extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
           Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Impulse response of  zp = "
             + String(dzp)));

    protected
      input Modelica_LinearSystems2.Types.TimeResponse response=
          Modelica_LinearSystems2.Types.TimeResponse.Impulse
        "type of time response";
      Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
          DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";

  algorithm
  // set sample time
      Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(
        dzp=dzp,
        tSpan=tSpan,
        response=response,
        x0=x0,
        defaultDiagram=defaultDiagram,
        device=device);

      annotation (interactive=true, Documentation(info="<html>


</html> "));
  end impulse;

  encapsulated function step
      "Step response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteZerosAndPoles dzp;
    input Real tSpan=0 "Simulation time span [s]";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Step response of  dzp = "
           + String(dzp)));

    protected
    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Step "type of time response";
    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";
  algorithm
   DiscreteZerosAndPoles.Plot.timeResponse(
      dzp=dzp,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

  equation

    annotation (interactive=true, Documentation(info="<html>

</html>"));
  end step;

  encapsulated function ramp
      "Ramp response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input DiscreteZerosAndPoles dzp;
    input Real tSpan=0 "Simulation time span [s]";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Ramp response of  dzp = "
           + String(dzp)));

    protected
    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Ramp "type of time response";

    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=zeros(
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp))
        "Initial state vector";
  algorithm
   Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(
      dzp=dzp,
      tSpan=tSpan,
      response=response,
      x0=x0,
      defaultDiagram=defaultDiagram,
      device=device);

    annotation (interactive=true, Documentation(info="<html>

</html> "));
  end ramp;

  encapsulated function initialResponse
      "Initial condition response plot of a discrete zeros-and-poles transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Types.TimeResponse;

      import Modelica_LinearSystems2.Utilities.Plot;

    input Modelica_LinearSystems2.DiscreteZerosAndPoles dzp;
    input Real tSpan=0 "Simulation time span [s]";

    input Modelica_LinearSystems2.Types.TimeResponse response=
        Modelica_LinearSystems2.Types.TimeResponse.Initial
        "type of time response";
    input Real y0 "Initial output (for initial condition plot)";

    extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
          Modelica_LinearSystems2.Internal.DefaultDiagramTimeResponse(heading="Initial response of  dzp = "
           + String(dzp) + "  with y0 = " + String(y0)));

    protected
    Modelica_LinearSystems2.DiscreteStateSpace dss=Modelica_LinearSystems2.DiscreteStateSpace(dzp);
    Real x0[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)]=
        Modelica.Math.Matrices.equalityLeastSquares(
        dss.A,
        fill(0, size(dss.B, 1)),
        dss.C,
        vector(y0)) "Initial state vector (for initial condition plot)";
  algorithm
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.timeResponse(
          dzp=dzp,
          tSpan=tSpan,
          response=response,
          x0=x0,
          defaultDiagram=defaultDiagram,
          device=device);

    annotation (interactive=true, Documentation(info="<html>


</html> "));
  end initialResponse;

  end Plot;

  encapsulated package Conversion
    import Modelica_LinearSystems2;
  function toDiscreteTransferFunction
      "Generate a DiscreteTransferFunction object from a DiscreteZerosAndPoles object"

      import Modelica;
      import Modelica_LinearSystems2.Math.Polynomial;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.Math.Complex;

    input DiscreteZerosAndPoles dzp
        "DiscreteZerosAndPoles transfer function of a system";
    output Modelica_LinearSystems2.DiscreteTransferFunction dtf;

    protected
    ZerosAndPoles zp=ZerosAndPoles(k=dzp.k, n1=dzp.n1, n2=dzp.n2, d1=dzp.d1, d2=dzp.d2);
    Real k;
    Complex z[:];
    Complex p[:];
    Polynomial pn;
    Polynomial pd;
  algorithm
    (z,p,k) := ZerosAndPoles.Analysis.zerosAndPoles(zp);
    pn := Polynomial(z)*Polynomial(k);
    pd := Polynomial(p);
    dtf.n := pn.c;
    dtf.d := pd.c;
    dtf.Ts := dzp.Ts;
    dtf.method := dzp.method;
    dtf.uName := dzp.uName;
    dtf.yName := dzp.yName;
    annotation (Documentation(info="<html>

</html>"));
  end toDiscreteTransferFunction;

  function toDiscreteStateSpace
      "Transform a DiscreteZerosAndPoles object into a DiscreteStateSpace object"
   //encapsulated function fromZerosAndPoles
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.Math.Vectors;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles.Internal;

    input DiscreteZerosAndPoles dzp
        "ZerosAndPoles transfer function of a system";
    output DiscreteStateSpace dss(
      redeclare Real A[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp),
        DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)],
      redeclare Real B[DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp),1],
      redeclare Real C[1,DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp)],
      redeclare Real D[1,1]) "Transfer function in StateSpace SISO form";

    protected
    Real A[2,2] "system matrix of partial 2nd order system";
    Real B[2,1] "input matrix of partial 2nd order system";
    Real C[1,2] "output matrix of partial 2nd order system";
    Real D[1,1] "feedthrough matrix of partial 2nd order system";
    Real a "system 'matrix' of partial 1st order system";
    Real b "input 'matrix' of partial 1st order system";
    Real c "output 'matrix' of partial 1st order system";
    Real d "feedthrough 'matrix' of partial 1st order system";
    Real aa "a2 + a1 +1";
    Real bb "b2 + b1 +1";
    Integer nx=max(DiscreteZerosAndPoles.Analysis.numeratorDegree(dzp),DiscreteZerosAndPoles.Analysis.denominatorDegree(dzp));
    Integer n_num1=size(dzp.n1, 1);
    Integer n_num2=size(dzp.n2, 1);
    Integer n_den1=size(dzp.d1, 1);
    Integer n_den2=size(dzp.d2, 1);
    Integer n_num=n_num1 + 2*n_num2;
    Integer n_den=n_den1 + 2*n_den2;

    Integer i_d=if n_num2 > n_den2 then 2*(n_num2 - n_den2) + 1 else 1;
    Integer i_k=if n_num2 > n_den2 then n_den2 - (n_num2 - n_den2) else n_den2;
    Integer i;
    Integer ili=if max(n_den2, n_num2) > 0 then i_d else  max(2, i_d);
    Real num[nx,2]=[dzp.n2; [dzp.n1,zeros(n_num1)]; zeros(max(0,nx - n_num2 - n_num1), 2)]
        "Numerator matrix, in order that indices are defined in all situations in all if clauses";
    Real den[nx,2]=[dzp.d2; [dzp.d1,zeros(n_den1)]; zeros(max(0,nx - n_den2 - n_den1), 2)]
        "Denominator matrix, in order that indices are defined in all situations in all if clauses";
    Real k[i_k + n_den1](each fixed=false)
        "Additional factors of the first and second order blocks, in order that the gain of the blocks is 1";
    Real k_total;
    Real eps=1e-5;//100*Modelica.Constants.eps;

    Boolean dZero=true;

  algorithm
    assert(n_num <= n_den,
      "DiscreteZerosAndPoles transfer function is not proper as required from StateSpace system:\n"
       + "  numerator degree (= " + String(n_num) +
      ") <= denominator degree (= " + String(n_den) + ") required.");

    if n_den > 0 then
      for i in 1:max(n_den2, n_num2) loop
        // State space systems of order 2
        if i <= n_den2 then
          if i <= n_num2 then
              // State space system in form (1)
            k[i] := Internal.scaleFactor2(
                num[i, 1],
                num[i, 2],
                den[i, 1],
                den[i, 2],eps);
          elseif 2*(i - n_num2) <= n_num1 then
              // State space system in form (1) with 2 first order numerator polynomials
            k[i] := Internal.scaleFactor2(
                num[2*(i - n_num2)-1, 1] + num[2*(i - n_num2), 1],
                num[2*(i - n_num2)-1, 1]*num[2*(i - n_num2), 1],
                den[i, 1],
                den[i, 2],eps);
          elseif  2*(i-n_num2) -1== n_num1 then
              // State space system in form (2) with 1 first order numerator polynomial
            k[i] := Internal.scaleFactor2(
                0,
                num[2*i-n_num2-1, 1],
                den[i, 1],
                den[i, 2],eps);
           else
              // State space system in form (3)
            k[i] := Internal.scaleFactor2(
                0,
                0,
                den[i, 1],
                den[i, 2],eps);
          end if;
        else
           // State space system in form (1) with 2 first order denominator polynomials

          k[i] := Internal.scaleFactor2(
              num[i, 1],
              num[i, 2],
              den[2*(i - n_den2)-1, 1] + den[2*(i - n_den2), 1],
              den[2*(i - n_den2)-1, 1]*den[2*(i - n_den2), 1],eps);
        end if;
      end for;

      for i in i_d:n_den1 loop
        // State space systems of order 1
        if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
           // State space system in form (4)
          k[i_k + i] := Internal.scaleFactor1(num[max(1, n_num2 + 2*(n_den2 -
            n_num2) + i), 1], den[n_den2 + i, 1],eps);
        elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
           // State space system in form (4)
          k[i_k + i] := Internal.scaleFactor1(num[max(1, n_num2 + i - i_d + 1),
            1], den[n_den2 + i, 1],eps);
        else
           // State space system in form (5)
          k[i_k + i] := Internal.scaleFactor1(0, den[n_den2 + i, 1],eps);
        end if;
      end for;

      k_total := dzp.k/product(k);

      dss.A := zeros(nx, nx);
      dss.B := zeros(nx, 1);
      dss.C := zeros(1, nx);
      dss.D := zeros(1, 1);

   // Calculation of matrices A, B, C, D
   //first elements of A, B, C and D

      if max(n_den2, n_num2) > 0 then

        A[1, :] := {0,1};
        B[1, 1] := 0;
            // Construct state space systems of order 2
        if 1 <= n_den2 then
          aa := num[1, 2] + num[1, 1] +1;
          bb := den[1, 2] + den[1, 1] +1;
          A[2, :] := {-den[1, 2],-den[1, 1]};
          B[2, 1] := if abs(bb)>eps and abs(aa)>eps then bb else 1;
          if 1 <= n_num2 then
                 // State space system in form (1)
            C := if abs(bb)>eps and abs(aa)>eps then k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]]/bb else k[1]*[num[1, 2] - den[1, 2],num[1, 1] - den[1, 1]];
            D := [k[1]];
            dZero := false;

          elseif 1 - n_num2 + 1 <= n_num1 then
                 // State space system in form (1) with 2 first order numerator polynomials
            aa := num[1, 1]*num[2, 1] + num[1, 1] + num[2, 1];
            B[2, 1] := if abs(bb)>eps  and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps and abs(aa)>eps then k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]]/bb else k[1]*[num[1, 1]*num[2, 1] - den[1, 2],num[1, 1] + num[2, 1] - den[1, 1]];
            D := [k[1]];
            dZero := false;
          elseif 1 - n_num2 == n_num1 then
                 // State space system in form (2) with 1 first order numerator polynomial
            aa := num[1, 1] + 1;
            B[2, 1] := if abs(bb)>eps  and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps  and abs(aa)>eps then k[1]*[num[1, 1],1]/bb else k[1]*[num[1, 1],1];
            D := [0];
            dZero := dZero and true;
          else
                 // State space system in form (3)
            B[2, 1] := if abs(bb)>eps  and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps  and abs(aa)>eps then  k[1]*[1,0]/bb else k[1]*[1,0];
            D := [0];
            dZero := dZero and true;
          end if;
       else
              // State space system in form (1) with 2 first order denominator polynomials
          bb := den[1, 1] + den[2, 1] + den[1, 1]*den[2, 1] + 1;
          A[2, :] := {-(den[1, 1]*den[2, 1]),-(den[1, 1] + den[2, 1])};
          B[2, 1] := if abs(bb)>eps then bb else 1;
          C := if abs(bb)>eps then  k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])]/bb else k[1]*[num[1, 2] - (den[1, 1]*den[2, 1]),num[1, 1] - (den[1, 1] + den[2, 1])];
          D := [k[1]];
          dZero := false;
        end if;

        dss.A[1:2, 1:2] := A;
        dss.B[1:2, 1] := vector(B);
        dss.C[1, 1:2] := vector(C);
        dss.D := D;

      else
        aa := num[1,1]+1;
        bb := den[1,1]+1;
        a := -den[1, 1];
        if 1 <= n_num1 then
          // State space system in form (4)
          b := if abs(bb)>eps then bb else num[1,1]-den[1,1];
          c := if abs(bb)>eps then k[1]*(num[1, 1] - den[1, 1])/bb else k[1];
          d := k[1];
          dZero := false;
        else
         // State space system in form (5)
          b := if abs(bb)>eps then bb else if n_num1>0 then num[1,1]-den[1,1] else 1;
          c := if abs(bb)>eps then k[1]/bb else k[1];
          d := 0;
          dZero := dZero and true;
        end if;
        dss.A[1, 1] := a;
        dss.B[1, 1] := b;
        dss.C[1, 1] := c;
        dss.D[1, 1] := d;

      end if;
    /// for i=2 to degree(system)
      A[1, :] := {0,1};
      B[1, 1] := 0;
      for i in 2:max(n_den2, n_num2) loop
           // Construct state space systems of order 2
        if i <= n_den2 then
          aa := num[i, 2] + num[i, 1] +1;
          bb := den[i, 2] + den[i, 1] +1;
          A[2, :] := {-den[i, 2],-den[i, 1]};
          B[2, 1] := if abs(bb)>eps and abs(aa)>eps then bb else 1;

          if i <= n_num2 then
                 // State space system in form (1)

            C := if abs(bb)>eps and abs(aa)>eps then k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]]/bb else k[i]*[num[i, 2] - den[i, 2],num[i, 1] - den[i, 1]];
            D := [k[i]];
            dZero := false;

          elseif 2*(i - n_num2) <= n_num1 then

            // State space system in form (1) with 2 first order numerator polynomials
            aa := num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] + num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] + 1;
            C := if abs(bb)>eps and abs(aa)>eps then k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]]/bb else k[i]*[num[2*i-n_num2-1, 1]*num[2*i-n_num2, 1] - den[i, 2],num[2*i-n_num2-1, 1] + num[2*i-n_num2, 1] - den[i, 1]];
            D := [k[i]];
            dZero := false;

          elseif 2*(i-n_num2) -1== n_num1 then
          // State space system in form (2) with 1 first order numerator polynomial
            aa := num[2*i-n_num2-1, 1]+1;
          B[2, 1] := if abs(bb)>eps and abs(aa)>eps then bb else 1;
            C := if abs(bb)>eps and abs(aa)>eps then k[i]*[num[2*i-n_num2-1, 1],1]/bb else k[i]*[num[2*i-n_num2-1, 1],1];
            D := [0];
            dZero := dZero and true;
          else
            // State space system in form (3)
            B[2, 1] := if abs(bb)>eps then bb else 1;
            C := if abs(bb)>eps then k[i]*[1,0]/bb else k[i]*[1,0];
            D := [0];
            dZero := dZero and true;
          end if;

        else
              // State space system in form (1) with 2 first order denominator polynomials
          bb := den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1] + den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1] + 1;
          A[2, :] := {-(den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),-(den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])};
          B[2, 1] := if abs(bb)>eps then bb else 1;
  //        C := if abs(bb)>eps then k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/den[max(2*(i-n_den2)-1,i), 1]/bb else k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])];
          C := if abs(bb)>eps then k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])]/bb else k[i]*[num[max(2*(i-n_num2),i), 2] - (den[max(2*(i-n_den2)-1,i), 1]*den[max(2*(i-n_den2),i), 1]),num[max(2*(i-n_num2),i), 1] - (den[max(2*(i-n_den2)-1,i), 1] + den[max(2*(i-n_den2),i), 1])];
          D := [k[i]];
          dZero := false;
        end if;
        dss.A[2*i, 1:2*i - 2] := B[2, 1]*dss.C[1, 1:2*i - 2];
        dss.A[2*i - 1:2*i, 2*i - 1:2*i] := A;
        dss.B[2*i, 1] := if dZero then 0 else B[2, 1]*dss.D[1, 1];
        dss.C[1, 1:2*i - 2] := if dZero then fill(0, 2*i - 2) else D[1, 1]*dss.C[
          1, 1:2*i - 2];
        dss.C[1, 2*i - 1:2*i] := vector(C);
        dss.D := D*dss.D;
      end for;
   //  for i in max(2,i_d):n_den1 loop
      for i in ili:n_den1 loop
           // Construct state space systems of order 1
        bb :=  den[n_den2 + i, 1] + 1;
        a := if abs(den[n_den2 + i, 1])>eps then -den[n_den2 + i, 1] else 0.0;

        if n_num2 <= n_den2 and 2*(n_den2 - n_num2) + i <= n_num1 then
              // State space system in form (4)

          c := if abs(bb)>eps then k[i_k + i]*(num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1] -  den[n_den2 + i, 1])/bb else 1.0;
  //        b := if abs(bb)>eps then bb else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1];
          b := if abs(bb)>eps then bb else num[max(1, n_num2 + 2*(n_den2 - n_num2) + i), 1]-  den[n_den2 + i, 1];
          d := k[i_k + i];
          dZero := false;
        elseif n_num2 > n_den2 and i - i_d + 1 <= n_num1 then
        // State space system in form (4)
          c := if abs(bb)>eps then k[i_k + i]*(num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1])/bb else 1.0;
  //        b := if abs(bb)>eps then bb else num[max(1, n_num2 + i - i_d + 1), 1];
          b := if abs(bb)>eps then bb else num[max(1, n_num2 + i - i_d + 1), 1] - den[n_den2 + i, 1];
          d := k[i_k + i];
          dZero := false;
        else
         // State space system in form (5)
          c := if abs(bb)>eps then k[i_k + i]/bb else k[i_k + i];
         b := if abs(bb)>eps then bb else 1;
  //        b := if abs(bb)>eps then den[n_den2 + i, 1] else 1;
          d := 0;
          dZero := dZero and true;
        end if;

  //###############################
  //        c := if abs(bb)>eps then k[1]*(num[1, 1] - den[1, 1])/bb else k[1];
  //        b := if abs(bb)>eps then bb else num[1,1]-den[1,1];
  //        d := k[1];
  //        dZero := false;
  //      else
  //       // State space system in form (5)
  //        c := if abs(bb)>eps then k[1]/bb else k[1];
  //        b := if abs(bb)>eps then bb else if n_num1>0 then num[1,1]-den[1,1] else 1;
  //        d := 0;
  //        dZero := dZero and true;
  //###############################

        dss.A[2*n_den2 + i, 1:2*n_den2 + i - 1] := b*dss.C[1, 1:2*n_den2 + i - 1];
        dss.A[2*n_den2 + i, 2*n_den2 + i] := a;
        dss.B[2*n_den2 + i, 1] := if dZero then 0 else b*dss.D[1, 1];
        dss.C[1, 1:2*n_den2 + i - 1] := if dZero then fill(0, 2*n_den2 + i - 1) else 
                d*dss.C[1, 1:2*n_den2 + i - 1];
        dss.C[1, 2*n_den2 + i] := c;
        dss.D := if dZero then [0] else d*dss.D;

        end for;

      dss.C := k_total*dss.C;
      dss.D := k_total*dss.D;
    else
      dss := DiscreteStateSpace(dzp.k);
    end if;

    dss.Ts := dzp.Ts;
    dss.method := dzp.method;

    annotation (overloadsConstructor=true, Documentation(info="<html>
<h4><font color=\"#008000\">Syntax</font></h4>
<table>
<tr> <td align=right>  ss </td><td align=center> =  </td>  <td> DiscreteZerosAndPoles.Conversion.toStateSpace<b>toStateSpace</b>(tf)  </td> </tr>
</table>
<h4><font color=\"#008000\">Description</font></h4>
<p>
This function transforms a zeros-poles-gain system representation into state space representation.
To achieve well numerical condition the DiscreteZerosAndPoles transfer function is transformed into state space
form by creating first and second order blocks that are connected
together in series. Every block is represented in controller
canonical form and scaled such that the gain from the input
of this block to its output is one (i.e. y(p=0) = u(p=0)),
if this is possible. Details are given below.
</p>
<b>Algorithmic details</b>
<p>
The DiscreteZerosAndPoles transfer function is defined as:
<blockquote><pre>
         product(p + n1[i]) * product(p^2 + n2[i,1]*p + n2[i,2])
  y = k*--------------------------------------------------------- * u
         product(p + d1[i]) * product(p^2 + d2[i,1]*p + d2[i,2])
</pre></blockquote>
<p>
This is treated as a series connection of first and second order
systems. If size(n1) == size(d1) and size(n2) == size(d2)
this gives the following sequence of operations:
</p>
<blockquote><pre>

        p^2 + n2[1,1]*p + n2[1,2]
  y_1 = ------------------------- * u
        p^2 + d2[1,1]*p + d2[1,2]
&nbsp;
        p^2 + n2[2,1]*p + n2[2,2]
  y_2 = ------------------------- * y_1
        p^2 + d2[2,1]*p + d2[2,2]
&nbsp;
     ...
&nbsp;
        p + n1[..]
  y_n = ---------- * y_(n-1)
        p + d1[..]
&nbsp;
    y = k*y_n

</pre></blockquote>

Based on this representation, evrey block with transfer function G(p) could be transformed into
<blockquote><pre>
  G(p) = k * F(p)

</pre>
with F(p) has unit gain. This leads to representations of the forms
</pre></blockquote>
<p>
<blockquote><pre>

           a2 + a1*p + p^2       a2      b2 + a1*b2/a2*p + b2/a2*p^2
  G(p) = -------------------- = ---- * ------------------------------ = k * F(p),  k = a2/b2  (1)

           b2 + b1*p + p^2       b2           b2 + b1*p + p^2
&nbsp;
for second order systems and
&nbsp;
           a + p     a     b + b/a*p
  G(p) = -------- = --- * ---------- = k * F(p),   k = a/b

           b + p     b      b + p
</p>
</pre></blockquote>
for first order systems respectively.

<p>
The complete system is now considered as the series connections of all the single unit gain transfer functions and an overall gain k with
<blockquote><pre>
  k = product(ki).

</pre></blockquote>
</p>


In the general case, the following system structures
and the corresponding state space systems can appear
(note, 'c' is the reciprocal local gain 1/k):
</p>
<blockquote><pre>
(1)
          a2 + a1*p + p^2           der(x1) = x2
    y = ---------------------  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
          b2 + b1*p + p^2                 y = c*((a2-b2)/b2*x1 + (a1-b1)/b2*x2 + u),  c = b2/a2
&nbsp;
(2)
             p + a                 der(x1) = x2
    y = ---------------- * u  -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                  y = k*(a/b2*x1 +x2/b2),  c = b2/a
&nbsp;
(3)
               1                  der(x1) = x2
    y = --------------- *u   -->  der(x2) = -b2*x1 - b1*x2 + b2*u
        b2 + b1*p + p^2                 y = c*x1/b2,  c = b2
&nbsp;
(4)
       a + p                       der(x) = -b*x + b*u
   y = ----- * u             -->        y = c*((a-b)/b*x + u),  c = b/a
       b + p
&nbsp;
(5)
         1
   y = ----- * u             -->   der(x) = -b*x + b*u
       b + p                            y = x,  c = b

</pre></blockquote>

If the sizes of the numerator and denominator polynomials
do not match, the small systems are built in the
following way:
</p>
<blockquote><pre>
(1) Build systems of form (1) by combining
    - 1 d2 and 1 n2
      (= 1 second order denominator and 1 second order numerator) or
    - 1 d2 and 2 n1 or
    - 2 d1 and 1 n2
(2) Build at most one system of form (2) by combining
    - 1 d2 and 1 n2
(3) Build systems of form (3) by
    - 1 d2
(4) Build systems of form (4) by combining
    - 1 d1 and 1 n1
(5) Build systems of form (5) by
    - 1 d1
</pre></blockquote>
<p>
The numeric properties of the resulting state space system
depends on which first and second order polynomials are
combined and connected together. From a numerical point of view, it
would therefore be useful to combine the polynomials
based on the numeric values of the polynomial coefficients,
(e.g., in a first step the polynomials could be sorted
according to their cut-off frequency).
</p>
<p>
However, this has the disadvantage that the structure of the
resulting state space system depends on the numeric
values of the polynomial coefficients. Since Modelica
environments perform symbolic pre-processing on equations,
this would mean that a change of a polynomial coefficient
requires to newly compile the state space system.
</p>
<p>
If, on the other hand, the structure of the state
space system depends only on dimension information
of the n1,n2,d1,d2 arrays, then the polynomial coefficients
can be changed without a new translation of the model.
This is the major reason why the structure of the
state space system in the implementation of this block
is based only on dimension information.
</p>
<p>
This is, e.g., not critical for the provided filters:
The dimension of the n1,n2,d1,d2 arrays depend for
filters only on the filter characteristics
(Bessel, Butterworth etc.), the filter type (low pass,
high pass etc.) and on the filter order. If any
of this data is changed, the model has to be
newly compiled. All the other filter data, such as
cut-off frequency or ripple amplitude, can be changed
without re-compilation of the model.
The DiscreteZerosAndPoles transfer function is now constructed
for the filters in such a way that the filter zeros
and poles are appropriately sorted to give better
numerical properties.
</p>
<p>
Another alternative implementation of the state
space system would be to use the function controller canonical
form that directly results from the transfer function.
The severe disadvantage
of this approach is that the structure of the state
space system from above is lost for the symbolic preprocessing.
If, e.g., index reduction has to be applied (e.g. since a
filter is used to realize a non-linear inverse model),
then the tool cannot perform the index reduction.
Example:
</p>
<p>
Assume, a generic first order state space system
is present
</p>
<blockquote><pre>
   <b>der</b>(x) = a*x + b*u
        y = c*x + d*u
</pre></blockquote>
<p>
and the values of the scalars a,b,c,d are parameters
that might be changed before the simulation starts.
If y has to be differentiated symbolically during code
generation, then
</p>
<blockquote><pre>
      <b>der</b>(y) = c*<b>der</b>(x) + d*<b>der</b>(u)
      <b>der</b>(x) = a*x + b*u
</pre></blockquote>
<p>
As a result, u needs to be differentiated too, and this
might not be possible and therefore translation might fail.
</p>
<p>
On the other hand, if the first order system is
defined to be a low pass filter and the state space
system is generated by keeping this structure, we have
(see form (5) above):
</p>
<blockquote><pre>
  <b>der</b>(x) = -b*x + u
        y = x
</pre></blockquote>
<p>
Differentiating y symbolically leads to:
</p>
<blockquote><pre>
     <b>der</b>(y) = <b>der</b>(x)
     <b>der</b>(x) = -b*x + u
</pre></blockquote>
<p>
Therefore, in this case, the derivative of u is not
needed and the tool can continue with the symbolic
processing.
</p>
 

</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   DiscreteZerosAndPoles q = Modelica_LinearSystems2.DiscreteZerosAndPoles.q();
   Modelica_LinearSystems2.DiscreteZerosAndPoles dzp=(q+1)/(q^2 + q + 1);

<b>algorithm</b>
  dss := Modelica_LinearSystems2.DiscreteZerosAndPoles.Conversion.toDiscreteStateSpace(dzp);
// dss.A = [0, 1; -1, -1],
// dss.B = [0; 3],
// dss.C = [0.333, 0.333],
// dss.D = [0],
// dss.B2 = [0; 0],

</pre></blockquote>

</html> "));
  end toDiscreteStateSpace;

  end Conversion;

  encapsulated package Import

  function fromModel
      "Generate a DiscreteZerosAndPoles record array from a state space representation resulted from linearization of a model"

      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteStateSpace;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;

      input String modelName "Name of the Modelica model";
      input Real T_linearize=0
        "point in time of simulation to linearize the model";
      input String fileName="dslin" "Name of the result file";
      input Modelica.SIunits.Time Ts=1 "Sample time";
      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";

    protected
      String fileName2=fileName + ".mat";
      Boolean OK1=simulateModel(problem=modelName, startTime=0, stopTime=T_linearize);
      Boolean OK2=importInitial("dsfinal.txt");
      Boolean OK3=linearizeModel(problem=modelName, resultFile=fileName, startTime=T_linearize, stopTime=T_linearize + 1);
      Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
      Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
      Integer nx=integer(nxMat[1, 1]);
      Integer nu=ABCDsizes[2] - nx;
      Integer ny=ABCDsizes[1] - nx;
      Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
      String xuyName[nx + nu + ny]=readStringMatrix(fileName2, "xuyName", nx + nu + ny);

      StateSpace ss(
        redeclare Real A[nx,nx],
        redeclare Real B[nx,nu],
        redeclare Real C[ny,nx],
        redeclare Real D[ny,nu]) "= model linearized at initial point";
      DiscreteStateSpace dss(
        redeclare Real A[nx,nx],
        redeclare Real B[nx,nu],
        redeclare Real C[ny,nx],
        redeclare Real D[ny,nu],
        redeclare Real B2[nx,nu]) "= model linearized at initial point";
      DiscreteStateSpace dss_siso(
        redeclare Real A[nx,nx],
        redeclare Real B[nx,1],
        redeclare Real C[1,nx],
        redeclare Real D[1,1],
        redeclare Real B2[nx,1]) "= model linearized at initial point";

        DiscreteZerosAndPoles dummy;

    public
      output DiscreteZerosAndPoles dzp[:,:];//=fill(dummy,ny,nu);
  algorithm
    ss.A := ABCD[1:nx, 1:nx];
    ss.B := ABCD[1:nx, nx + 1:nx + nu];
    ss.C := ABCD[nx + 1:nx + ny, 1:nx];
    ss.D := ABCD[nx + 1:nx + ny, nx + 1:nx + nu];
    ss.uNames := xuyName[nx + 1:nx + nu];
    ss.yNames := xuyName[nx + nu + 1:nx + nu + ny];
    ss.xNames := xuyName[1:nx];

    dss := DiscreteStateSpace(ss, Ts=Ts, method=method);

    dzp := DiscreteStateSpace.Conversion.toDiscreteZerosAndPolesMIMO(dss);

  //   for ic in 1:ny loop
  //     for ib in 1:nu loop
  //       dss_siso := DiscreteStateSpace(
  //         A=dss.A,
  //         B=matrix(dss.B[:, ib]),
  //         C=transpose(matrix(dss.C[ic, :])),
  //         D=matrix(dss.D[ic, ib]),
  //         B2=matrix(dss.B2[:, ib]),
  //         Ts=dss.Ts,
  //         method=dss.method);
  //       dzp[ic, ib] := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(
  //         dss_siso);
  //     end for;
  //   end for;

  //
  // //    zp := StateSpace.Conversion.toZerosAndPolesMIMO(result);
  //     for i in 1:ny loop
  //       for j in 1:nu loop
  //         dzp[i,j] := DiscreteZerosAndPoles(zp=zp[i,j], Ts=Ts, method=method);
  //       end for;
  //     end for;

      annotation (interactive=true, Documentation(info="function fromModel 
  \"Generate a ZerosAndPoles record array from a state space representation resulted from linearization of a model\"

  import Modelica;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;

  input String modelName \"Name of the Modelica model\" annotation(Dialog(translatedModel));
  input Real T_linearize=0 \"point in time of simulation to linearize the model\";
  input String fileName=\"dslin\" \"Name of the result file\";

  annotation (interactive=true, Documentation(info=\"<html>

 
</html> \"));

"));
  end fromModel;

    encapsulated function fromFile
      "Generate a DiscreteZerosAndPoles record by reading the polynomial coefficients or zeros and poles from a file"
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica;
      import Modelica_LinearSystems2.DataDir;

      input String fileName=DataDir + "/dzp.mat"
        "Name of the discrete zeros and poles data file"        annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                        caption="state space system data file")));

    protected
      input Integer n1n2d1d2[4]=if ZerosAndPoles.Internal.checkRepresentation(
          fileName) then ZerosAndPoles.Internal.numberOfRealZerosAndPoles_zp(
          fileName) else ZerosAndPoles.Internal.numberOfRealZerosAndPoles_pc(
          fileName);
      input Integer n1=n1n2d1d2[1];
      input Integer n2=n1n2d1d2[2];
      input Integer d1=n1n2d1d2[3];
      input Integer d2=n1n2d1d2[4];
      input Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2];
      input Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
      output DiscreteZerosAndPoles dzp(n1=fill(0, n1), n2=fill(0, n2, 2), d1=fill(0, d1), d2=fill(0, d2, 2));
    algorithm
    //Whenever this function becomes operational the code must be rewritten if fromFile_pc2 and fromFile_zp2 are in the 'constructor'

      dzp := if ZerosAndPoles.Internal.checkRepresentation(fileName) then DiscreteZerosAndPoles.Internal.fromFile_zp( fileName) else DiscreteZerosAndPoles.Internal.fromFile_pc(
        fileName);

      annotation (Documentation(info="<html>

</html> "));
    end fromFile;

  end Import;

  encapsulated package Internal
    "Internal library of record Filter (should not be directly used by user)"

    import Modelica;
    import Modelica_LinearSystems2;
    extends Modelica.Icons.Library;

    function numberOfRealZeros2 "Calculate number of real zeros"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

    input Modelica_LinearSystems2.DiscreteTransferFunction dtf
        "DiscreteTransferFunction";
    output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(dtf.n)));
    algorithm
    end numberOfRealZeros2;

    function numberOfRealPoles "Calculate number of real poles"
      import Modelica;
      import Modelica_LinearSystems2.Internal;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

    input Modelica_LinearSystems2.DiscreteTransferFunction dtf
        "TransferFunction";
      output Integer result=Internal.numberOfRealZeros(Polynomial.roots(Polynomial(
          dtf.d)));
    algorithm
    end numberOfRealPoles;

  encapsulated function scaleFactor1
      "Return scale factor for first order block"
      import Modelica;
    input Real n "(z+n)/(z+d)";
    input Real d "(z+n)/(z+d)";
    input Real small=100*Modelica.Constants.eps;
    output Real k
        "= (d+1)/(n+1), if d+1, n+1 are not zero, otherwise special cases";
  algorithm
      k := if abs(d+1) > small  and abs(n+1) > small then abs(d+1)/abs(n+1) else 1;
  end scaleFactor1;

    function scaleFactor2 "Return scale factor for second order block"
      import Modelica;
    input Real n1 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real n2 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real d1 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real d2 "(z^2 + n1*z + n2)/(z^2 + d1*z + d2)";
    input Real small=100*Modelica.Constants.eps;
    output Real k
        "= (d2+d1+1)/(n2+n1+1), if numerator and denominator are not zero, otherwise special cases";
    algorithm
      k := if abs(d2+d1+1) > small and abs(n2+n1+1) > small then (d2+d1+1)/(n2+n1+1) else 1;
    end scaleFactor2;

  function findMatrixName
    input String filename;
    input String matName="z";

    output Integer m;

  external "C" m=  findMatrixName(
          filename,
          matName,
          "NoClass");

    annotation (Include="#include <matrixop.h> 
#include <matrixop.c> 
 
 
#if !defined(DYMOLA_DSPACE) && !defined(NO_FILE)
#include <amat.h>
#include <sprwat.h>
#endif
 
extern int findMatrixName(const char* fil,const char* matname, char *noClass) {
int found=0;
 
 
#if !defined(DYMOLA_DSPACE) && !defined(NO_FILE)
{
        AmatGetFile afile;
        Amatrix amatrix;
 
        int ret=amatGetOpen((char*)fil,noClass,(char*)0,&afile);
 
        Assert(ret==0,amatError);
        for(;ret==0 && !found;)  
{
                amatInit(&amatrix);
                ret=amatGetMatrix(&afile, &amatrix);
                if (ret<=1 && strcmp(matname,amatrix.name)==0)
                  found=1;
                else
                  found=0;
 
                amatDel(&amatrix);
        }
        amatGetClose(&afile);
}
#else
        Assert(false, 'nn');
#endif
        return found;
}");
  end findMatrixName;

  function fromRealAndImag
      "Generate a complex vector from a real part vector and imaginary part vector "

      import Modelica_LinearSystems2.Math.Complex;

    input Real real[:];
    input Real imag[size(real, 1)];
    output Complex result[size(real, 1)] "Number of real zeros";
  algorithm
    for i in 1:size(real, 1) loop
      result[i].re := real[i];
      result[i].im := imag[i];
    end for;
  end fromRealAndImag;

  encapsulated function fromFile_pc
      "Generate a zeros and poles data record by reading the polynomial coefficients from a file (default file name is pc.mat)"
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Math.Complex;
      import Modelica;

    input String fileName="pc.mat" "Name of the zeros and poles data file" 
                                                     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));

    protected
    input Integer n1n2d1d2[4]=ZerosAndPoles.Internal.numberOfRealZerosAndPoles_pc(fileName);
    input Integer n1=n1n2d1d2[1];
    input Integer n2=n1n2d1d2[2];
    input Integer d1=n1n2d1d2[3];
    input Integer d2=n1n2d1d2[4];
    input Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2];
    input Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
    output DiscreteZerosAndPoles dzp(
      n1=fill(0, n1),
      n2=fill(
            0,
            n2,
            2),
      d1=fill(0, d1),
      d2=fill(
            0,
            d2,
            2));

    protected
    Integer n1_2=if n1 > 0 then 1 else 0 "second dimension of n1-matrix";
    Integer n2_2=if n2 > 0 then 2 else 0 "second dimension of n2-matrix";
    Integer d1_2=if d1 > 0 then 1 else 0 "second dimension of d1-matrix";
    Integer d2_2=if d2 > 0 then 2 else 0 "second dimension of d2-matrix";

    Real k=scalar(readMatrix(
          fileName,
          "k",
          1,
          1));
    Real n1Vector[n1]=vector(readMatrix(
          fileName,
          "n1",
          n1,
          n1_2)) "coefficients of first order numenator polynomials";
    Real n2Matrix[n2,n2_2]=readMatrix(
          fileName,
          "n2",
          n2,
          n2_2) "coefficients of second order denominator polynomials";
    Real d1Vector[d1]=vector(readMatrix(
          fileName,
          "d1",
          d2,
          d1_2)) "coefficients of first order denominator polynomials";
    Real d2Matrix[d2,d2_2]=readMatrix(
          fileName,
          "d2",
          d2,
          d2_2) "coefficients of second order numenator polynomials";
    Real Ts[1,1]=readMatrix(fileName, "Ts", 1, 1);

  algorithm
    dzp.k := k;
    dzp.n1 := if n1 > 0 then n1Vector else fill(0, 0);
    dzp.n2 := if n2 > 0 then n2Matrix else fill(0, 0, 2);
    dzp.d1 := if d1 > 0 then d1Vector else fill(0, 0);
    dzp.d2 := if d2 > 0 then d2Matrix else fill(0, 0, 2);
    dzp.Ts := scalar(Ts);

  end fromFile_pc;

  encapsulated function fromFile_zp
      "Generate a zeros and poles data record by reading poles and zeros from a file (default file name is zp.mat)"

      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.ZerosAndPoles;
      import Modelica_LinearSystems2.DiscreteZerosAndPoles;
      import Modelica_LinearSystems2.Math.Complex;

    input String fileName="zp.mat" "Name of the zeros and poles data file" 
                                                     annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
                      caption="state space system data file")));
    protected
    input Integer n1n2d1d2[4]=
        ZerosAndPoles.Internal.numberOfRealZerosAndPoles_zp(fileName);
    input Integer n1=n1n2d1d2[1];
    input Integer n2=n1n2d1d2[2];
    input Integer d1=n1n2d1d2[3];
    input Integer d2=n1n2d1d2[4];
    input Integer zSize=n1n2d1d2[1] + 2*n1n2d1d2[2];
    input Integer pSize=n1n2d1d2[3] + 2*n1n2d1d2[4];
    public
    output DiscreteZerosAndPoles dzp(
      n1=fill(0, n1),
      n2=fill(
            0,
            n2,
            2),
      d1=fill(0, d1),
      d2=fill(
            0,
            d2,
            2));

    protected
    Integer z_2=if zSize > 0 then 2 else 0 "second dimension of zeros-matrix";
    Integer p_2=if pSize > 0 then 2 else 0 "second dimension of poles-matrix";

    Real k=scalar(readMatrix(
          fileName,
          "k",
          1,
          1));
    Real zerosMatrix[zSize,z_2]=readMatrix(
          fileName,
          "z",
          zSize,
          z_2) "zeros in rows of real parts and imaginary parts";
    Real polesMatrix[pSize,p_2]=readMatrix(
          fileName,
          "p",
          pSize,
          p_2) "poles in rows of real parts and imaginary parts";
    Real Ts[1,1]=readMatrix(fileName, "Ts", 1, 1);

    Complex zeros[:]=if zSize > 0 then ZerosAndPoles.Internal.fromRealAndImag(
        zerosMatrix[:, 1], zerosMatrix[:, z_2]) else fill(Complex(0), 0);
    Complex poles[:]=if pSize > 0 then ZerosAndPoles.Internal.fromRealAndImag(
        polesMatrix[:, 1], polesMatrix[:, p_2]) else fill(Complex(0), 0);

  algorithm
    dzp := DiscreteZerosAndPoles(
        k=k,
        z=zeros,
        p=poles,
        Ts=scalar(Ts));
  end fromFile_zp;

  end Internal;

  annotation (
    defaultComponentName="filter",
    Documentation(info="<html>
<p>
This record defines a discrete transfer function by its zeros, poles and a gain:
</p>
<pre>         product(q - z[i])
  y = k*------------------- * u
         product(q - n[i])
</pre>
<p>
where z[:] is a Complex vector of zeros, n[:] is a Complex
vector of poles and k is an additional multiplicative factor.
The elements of the two Complex vectors must either be real
numbers or conjugate complex pairs (in order that their product
results in a polynomial with Real coefficients).
</p>
<p>
In the record, the zeros and poles are transformed
into a product of first and second order polynomials.

</p>
<p>
A DiscreteZeroAndPole transfer function is internally stored by the coefficients
of first and second order polynomials, and by an additional
multiplicative factor k:
</p>
<pre>         product(q + n1[i]) * product(q^2 + n2[i,1]*q + n2[i,2])
  y = k*---------------------------------------------------------
         product(q + d1[i]) * product(q^2 + d2[i,1]*q + d2[i,2])
</pre>
<p>
Note, the degrees of the numerator and denominator
polynomials are given as:
</p>
<pre>
   degree of numerator   = size(n1,1) + 2*size(n2,1);
   degree of denominator = size(d1,1) + 2*size(d2,1);
</pre>
<p>
Example:
</p>
<pre>                          (q+0.5)
  dzp = 4* -------------------------------------
           (p - 0.5)*(p - (0.6+j*0.3))*(p - (0.6-j*0.3))
</pre>
<p>
with j = Complex.j(); is defined as
</p>
<pre> 
   <b>import</b> Modelica_LinearSystems2.Math.Complex; 
   <b>import</b> Modelica_LinearSystems2.DiscreteZerosAndPoles;
   j = Complex.j();
   
   dzp = ZerosAndPoles(z = {Complex(-0.5)},
                       p = {Complex(0.5),
                            0.6+j*0.3,
                            0.6-j*0.3}, 
                            k=4);
</pre>
</html>"),
    DymolaStoredErrors);
end DiscreteZerosAndPoles;
