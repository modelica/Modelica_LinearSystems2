within Modelica_LinearSystems2.Utilities.Plot;
function parameterizedCurves
  "Plot parametrized curve with one or more branches"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Utilities.Plot.Internal;

  input Modelica_LinearSystems2.Utilities.Plot.Records.ParametrizedCurves diagram
    "Parametrized curve data points" annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
    Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);

protected
  Boolean ok "True, if call is ok";
  Real mmToPixel= device.windowResolution/25.4;
  Real position[4];
  Integer one=1;
  Integer k;
  Integer j;
  Integer id;
  Integer nProperties=size(diagram.curveProperties,1);
  Integer nBranches=size(diagram.X,1);
  Integer colors[nBranches,3] "Line colors";
  Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern patterns[nBranches]
    "Line patterns";
  Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol symbols[nBranches]
    "Line symbols";
  Real thicknesses[nBranches] "Line thicknesses";
algorithm
  // Create diagram
  position := {device.xTopLeft,
               device.yTopLeft,
               device.diagramWidth,
               diagram.heightRatio*device.diagramWidth}*mmToPixel;

  id:= DymolaCommands.Plot.createPlot(
    id=-1,
    position=integer(position),
    erase=true,
    autoscale=true,
    autoerase=false,
    subPlot=1,
    heading=diagram.heading,
    grid=diagram.grid,
    logX=diagram.logX,
    logY=diagram.logY,
    bottomTitle=diagram.xLabel,
    leftTitle=diagram.yLabel,
    color=device.autoLineColor,
    legend=diagram.legend,
    legendHorizontal=diagram.legendHorizontal,
    legendFrame=diagram.legendFrame,
    legendLocation=Integer(diagram.legendLocation));

  // Plot parameterized curve
  if nProperties == 0 then
    DymolaCommands.Plot.plotParametricCurves(
      x=diagram.X,
      y=diagram.Y,
      s=diagram.s,
      xName=diagram.xName,
      yName=diagram.yName,
      sName=diagram.sName,
      legends=diagram.legends,
      id = id,
      labelWithS=diagram.labelWithS);
  else
    for i in 1:nBranches loop
      k := i;
      j :=mod(k, nProperties) + 1
        "if k is replaced by i, Dymola gives an error about assignment of Real to Integer";
      colors[i,:]    :=diagram.curveProperties[j].lineColor;
      patterns[i]    :=diagram.curveProperties[j].linePattern;
      symbols[i]     :=diagram.curveProperties[j].lineSymbol;
      thicknesses[i] :=diagram.curveProperties[j].lineThickness;
    end for;
    DymolaCommands.Plot.plotParametricCurves(
      x=diagram.X,
      y=diagram.Y,
      s=diagram.s,
      xName=diagram.xName,
      yName=diagram.yName,
      sName=diagram.sName,
      legends=diagram.legends,
      id=id,
      labelWithS=diagram.labelWithS,
      colors=colors,
      patterns=Internal.convertToDymolaPattern(patterns),
      markers=Internal.convertToDymolaMarker(symbols),
      thicknesses=thicknesses);
  end if;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Plot.<strong>parameterizedCurves</strong>(diagram, device)
</pre></blockquote>

<h4>Description</h4>
<p>
This function plots a set of parametrized curves that depend on the same path parameter s
in one window. The set of parameterized curves is defined by:
</p>
<pre>
   s = s[i]   // Vector of s-values
   X = X[j,i] // X=X(s), where s[i] is the s-value and j is the s-branch of the x-value
   Y = Y[j,i] // Y=Y(s), where s[i] is the s-value and j is the s-branch of the y-value
</pre>

<h4>Example</h4>
<p>
Execute function <a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.plotParameterizedCurve1\">Examples.plotParameterizedCurve1</a>
that is defined as parameterized sine and cosine-functions:
</p>
<blockquote><pre>
s = linspace(0, Modelica.Units.Conversions.from_deg(300), 100);
<strong>for</strong> i <strong>in</strong> 1:nPoints <strong>loop</strong>
  X[1, i] := cos(s[i]);
  Y[1, i] := sin(s[i]);
  X[2, i] := X[1, i] - 0.5;
  Y[2, i] := Y[1, i];
<strong>end for</strong>;
Plot.parameterizedCurves(diagram=
  Plot.Records.ParametrizedCurves(
    X=X,
    Y=Y,
    s=s));
</pre></blockquote>
<p>
This yields the following diagram (the menu on the right lower part is displayed when moving
the cursor on one curve point; then all points belonging to the same path parameter value are
marked with a red square):
</p>
<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotParameterizedCurve1.png\"/>
</div>
</html>"));
end parameterizedCurves;
