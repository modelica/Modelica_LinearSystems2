within Modelica_LinearSystems2.Utilities.Plot.Examples;
function plotParameterizedCurve2 "Plot two circles as parameterized curve"
  extends Modelica.Icons.Function;

  import Modelica.SIunits.Conversions.from_deg;
  import Modelica.Math.sin;
  input Modelica.SIunits.Conversions.NonSIunits.Angle_deg maxAngle = 300
    "Maximum opening angle";
  input Integer nPoints(min=2) = 100 "Number of points";
protected
  Real s[nPoints] = linspace(0,from_deg(maxAngle),nPoints);
  Real X[2,nPoints];
  Real Y[2,nPoints];
algorithm
  for i in 1:nPoints loop
    X[1,i] :=cos(s[i]);
    Y[1,i] :=sin(s[i]);
    X[2,i] :=X[1,i] - 0.5;
    Y[2,i] :=Y[1,i];
  end for;
  Modelica_LinearSystems2.Utilities.Plot.parameterizedCurves(
    diagram = Modelica_LinearSystems2.Utilities.Plot.Records.ParametrizedCurves(
      X=X,
      Y=Y,
      s=s,
      xName="x",
      heading="Two circles as parameterized curves",
      legends={"circle1","circle2"},
      labelWithS=true,
      curveProperties={
        Modelica_LinearSystems2.Utilities.Plot.Records.CurveProperties(
          lineColor={0,127,0},
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.DashDot,
          lineThickness=0.5),
        Modelica_LinearSystems2.Utilities.Plot.Records.CurveProperties(
          lineColor={170,85,255},
          linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid,
          lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Cross)},
      xLabel="x(s)",
      yLabel="y(s)",
      legend=true));

  annotation (
    Documentation(info="<html>
<p>
This function plots the following parameterized cuves diagram using function
<a href=\"Modelica_LinearSystems2.Utilities.Plot.parameterizedCurves\">parameterizedCurves</a>
</p>

<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotParameterizedCurve2.png\"/></p>
</html>"));
end plotParameterizedCurve2;
