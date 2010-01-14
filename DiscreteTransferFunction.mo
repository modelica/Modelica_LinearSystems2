within Modelica_LinearSystems2;
record DiscreteTransferFunction
  "Discrete transfer function description of a single input, single output system (data + operations)"

  extends Modelica.Icons.Record;
  import Modelica_LinearSystems2.Math.Polynomial;
  import Modelica_LinearSystems2;

  Real n[:] "Coefficients of numerator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));
  Real d[:] "Coefficients of denominator polynomial (in descending order)" annotation(Dialog(group="y = n*{s^m, ... , s, 1} / (d*{s^r, ... , s, 1}) * u"));

  Modelica.SIunits.Time Ts "Sample time" 
       annotation(Dialog(group="Data used to construct discrete from continuous system"));

  Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
    "Discretization method" 
        annotation(Dialog(group="Data used to construct discrete from continuous system"));

  String uName="u" "Name of input signal"    annotation(Dialog(group="Signal names"));
  String yName="y" "Name of output signal"  annotation(Dialog(group="Signal names"));

/* If the numerator polynomial has no coefficients, the transfer function
   is zero. The denominator polynomial must always have at
   least one coefficient, such as {1}
*/

  encapsulated operator 'constructor'
    "Default constructor for a transfer function"
    import Modelica;
    import Modelica_LinearSystems2.TransferFunction;

    function fromReal
      "Generate a DiscreteTransferFunction data record from a Real value"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteTransferFunction;

      input Real r "Value of Real variable";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(n={r}, d={1});

    algorithm
      dtf.uName := uName;
      dtf.yName := yName;
      annotation (overloadsConstructor=true);
    end fromReal;

      encapsulated function fromArrays
      "Generate a DiscreteTransferFunction data record from numerator and denominator array"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2;

           input Real n[:] "Coefficients of numerator polynomial";
           input Real d[:] "Coefficients of denominator polynomial";
           input Modelica.SIunits.Time Ts "Sample time";
           input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method";

           input String uName = "" "input name";
           input String yName = "" "output name";

           output DiscreteTransferFunction dtf(redeclare Real n[size(n, 1)], redeclare
          Real d[                                                                             size(d, 1)])
        "Transfer function";

      algorithm
                 //this is the constructor algorithm
           assert(size(d, 1) > 0, "Input denominator d must have at least one element, however\n"
              + "d is an empty vector");
           dtf.n := n;
           dtf.d := d;
           dtf.Ts := Ts;
           dtf.method := method;
           dtf.uName := uName;
           dtf.yName := yName;
      end fromArrays;

    function fromPolynomials
      "Generate a DiscreteTransferFunction data record from a numerator and denominator polynomial"
      import Modelica;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.Math.Polynomial;

      input Polynomial n "Numerator polynomial";
      input Polynomial d "Denominator polynomial";
      input String uName="" "input name";
      input String yName="" "output name";
      output DiscreteTransferFunction dtf(n=n.c, d=d.c,uName=uName, yName=yName);

    algorithm
      annotation (overloadsConstructor=true);
    end fromPolynomials;

    function fromTransferFunction
      "Generate a DiscreteTransferFunction data record from a continuous Transfer function"
      import Modelica;
      import Modelica_LinearSystems2;
      import Modelica_LinearSystems2.DiscreteTransferFunction;
      import Modelica_LinearSystems2.TransferFunction;
      import Modelica_LinearSystems2.StateSpace;
      import Modelica_LinearSystems2.DiscreteStateSpace;

      input TransferFunction tf "continuous transfer function";
      input Modelica.SIunits.Time Ts "Sample time" 
           annotation(Dialog(group="Data used to construct discrete from continuous system"));

      input Modelica_LinearSystems2.Types.Method method=Modelica_LinearSystems2.Types.Method.Trapezoidal
        "Discretization method" 
            annotation(Dialog(group="Data used to construct discrete from continuous system"));

      output DiscreteTransferFunction dtf;
    protected
      StateSpace ss = StateSpace(tf);
      DiscreteStateSpace dss = DiscreteStateSpace(ss,Ts,method);

    algorithm
      dtf := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss);
      annotation (overloadsConstructor=true);
    end fromTransferFunction;
  end 'constructor';

    encapsulated operator function 'String'
    "Transform TransferFunction into a String representation"
    import Modelica_LinearSystems2.Math.Polynomial;
    import Modelica_LinearSystems2.DiscreteTransferFunction;

    input DiscreteTransferFunction dtf
      "Discrete transfer function to be transformed in a String representation";
    input Integer significantDigits=6
      "Number of significant digits that are shown";
    input String name="z" "Independent variable name used for printing";
    output String z="";
  protected
    Integer n_num=size(dtf.n, 1) - 1;
    Integer n_den=size(dtf.d, 1) - 1;
    Boolean numParenthesis;
    algorithm
    if n_num == -1 then
      z := "0";
    else
      numParenthesis := n_num > 0 and not (n_den == 0 and dtf.d[1] == 1);
      if numParenthesis then
        z := "(";
      end if;
       z := z + String(
            Polynomial(dtf.n),
            significantDigits,
            name);

      if numParenthesis then
        z := z + ")";
      end if;
    if n_den > 0 or dtf.d[1] <> 1 then
        if n_den > 0 then
          z := z + "/(";
        else
          z := z + "/";
        end if;

        z := z + String(
              Polynomial(dtf.d),
              significantDigits,
              name);

        if n_den > 0 then
          z := z + ")";
        end if;
      end if;
    end if;
    end 'String';

encapsulated package Plot "Functions to plot state space system responses"

encapsulated function bode "Plot transfer function as bode plot"
  import Modelica;
  import Modelica.Utilities.Strings;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Internal;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.DiscreteTransferFunction;
  import Modelica_LinearSystems2.Math.Complex;

  import Modelica_LinearSystems2.Utilities.Plot;

  import SI = Modelica.SIunits;

  input DiscreteTransferFunction dtf "DiscreteTransfer function to be plotted";
  input Integer nPoints(min=2) = 200 "Number of points";
  input Boolean autoRange=true
        "= true, if abszissa range is automatically determined";
  input Modelica.SIunits.Frequency f_min(min=0) = 0.1
        "Minimum frequency value, if autoRange = false" 
                                                    annotation(Dialog(enable=not autoRange));
  input Modelica.SIunits.Frequency f_max(min=0) = 10
        "Maximum frequency value, if autoRange = false"                                              annotation(Dialog(enable=not autoRange));

  input Boolean magnitude=true "= true, to plot the magnitude of tf" 
                                                                    annotation(choices(__Dymola_checkBox=true));
  input Boolean phase=true "= true, to plot the pase of tf" annotation(choices(__Dymola_checkBox=true));

  extends Modelica_LinearSystems2.Internal.PartialPlotFunction(defaultDiagram=
        Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot(heading="Bode plot of  dtf = "
         + String(dtf)));

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
  TransferFunction tf=TransferFunction(n=dtf.n, d=dtf.d);

  Plot.Records.Curve curves[2];
  Integer i;
  Plot.Records.Diagram diagram2[2];

algorithm
        // Determine frequency vector f
  if autoRange then
    (numZerosZ,denZerosZ) := TransferFunction.Analysis.zerosAndPoles(tf);
  else
    numZerosZ := fill(Complex(0), 0);
    denZerosZ := fill(Complex(0), 0);
  end if;

  numZeros := fill(Complex(0),0);
 // numZeros := fill(Complex(0),size(numZerosZ,1));
  // for i in 1:size(numZerosZ,1) loop
  //   numZeros[i] := Complex.log(numZerosZ[i])/dtf.Ts;
  // end for;

  denZeros := fill(Complex(0),size(denZerosZ,1));
  for i in 1:size(denZerosZ,1) loop
    denZeros[i] := Complex.log(denZerosZ[i])/dtf.Ts;
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
    z[i] := Complex.exp(Complex(0,w[i]*dtf.Ts));
    c := TransferFunction.Analysis.evaluate(
          tf,
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
<h4><font color=\"#008000\">Syntax</font></h4>
<blockquote><pre>
TransferFunction.Plot.<b>plotBode</b>(tf)
   or
TransferFunction.Plot.<b>plotBode</b>(tf, nPoints, autoRange, f_min, f_max, magnitude=true, phase=true, defaultDiagram=<a href=\"Modelica://Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot\">Modelica_LinearSystems2.Internal.DefaultDiagramBodePlot</a>(), device=<a href=\"Modelica://Modelica_LinearSystems2.Utilities.Plot.Records.Device\">Modelica_LinearSystems2.Utilities.Plot.Records.Device</a>() )
</pre></blockquote>
<h4><font color=\"#008000\">Description</font></h4>
<p>
Plots the bode-diagram of a transfer function.


</p>

<h4><font color=\"#008000\">Example</font></h4>
<blockquote><pre>
   TransferFunction s = Modelica_LinearSystems2.TransferFunction.s();  
   Modelica_LinearSystems2.TransferFunction tf =(s^2 + 5*s + 7)/(s^2 + 5*s + 6);
   
<b>algorithm</b>
   Modelica_LinearSystems2.TransferFunction.Plot.plotBode(tf)
//  gives:
</pre></blockquote>

</p>
 
<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodeMagnitude.png\">
<br>
<img src=\"modelica://Modelica_LinearSystems2/Extras/Images/bodePhase.png\">
 
</blockquote>
<p>


</html> "));
end bode;

end Plot;

encapsulated package Analysis

end Analysis;

end DiscreteTransferFunction;
