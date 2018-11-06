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

  id:= DymolaCommands.Plot.createPlot(id=-1,
                  position=integer(position),
                  y={""},
                  heading=diagram.heading,
                  erase=true,
                  autoscale=true,
                  autoerase=false,
                  grid=diagram.grid,
                  legend=diagram.legend,
                  legendLocation=Integer(diagram.legendLocation),
                  legendHorizontal=diagram.legendHorizontal,
                  legendFrame=diagram.legendFrame,
                  logX=diagram.logX,
                  logY=diagram.logY,
                  legends={""},
                  subPlot=1,
                  leftTitle=diagram.yLabel,
                  bottomTitle=diagram.xLabel,
                  colors=fill({-1, -1, -1}, 1),
                  patterns={LinePattern.Solid},
                  markers={MarkerStyle.None},
                  thicknesses={0.5},
                  axes={1},
                  displayUnits={""});

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
                  labelWithS=diagram.labelWithS,
                  colors=fill({-1, -1, -1}, size(diagram.Y, 1)),
                  patterns=fill(LinePattern.Solid, size(diagram.Y, 1)),
                  markers=fill(MarkerStyle.None, size(diagram.Y, 1)),
                  thicknesses=fill(0.25, size(diagram.Y, 1)),
                  axes=fill(1, size(diagram.Y, 1)));
  else
    for i in 1:nBranches loop
      k := i;
      j := mod(k, nProperties) + 1
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
      thicknesses=thicknesses,
      axes=fill(1, size(diagram.Y, 1)));
  end if;

/*
function plotParametricCurves "plot parametric curves"
  input Real x[:, size(s, 1)] "x(s) vectors";
  input Real y[size(x, 1), size(s, 1)] "y(s) vectors";
  input Real s[:] "s values";
  input String xName := "" "The name of the x variable";
  input String yName := "" "The name of the y variable";
  input String sName := "" "The name of the s parameter";
  input String legends[:] "Legends describing plotted data";
  input Integer id := 0 "Identity of window (0-means last)";
  input Integer colors[size(y, 1), 3] "Line colors";
  input Integer patterns[size(y, 1)] "Line patterns, e.g., LinePattern.Solid";
  input Integer markers[size(y, 1)] "Line markers, e.g., MarkerStyle.Cross";
  input Real thicknesses[size(y, 1)] "Line thicknesses";
  input Boolean labelWithS := false "if true, output values of s along the curve";
*/

  annotation (
    Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Plot.<b>parameterizedCurves</b>(diagram, device)
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
s = linspace(0, Modelica.SIunits.Conversions.from_deg(300), 100);
<b>for</b> i <b>in</b> 1:nPoints <b>loop</b>
  X[1, i] := cos(s[i]);
  Y[1, i] := sin(s[i]);
  X[2, i] := X[1, i] - 0.5;
  Y[2, i] := Y[1, i];
<b>end for</b>;
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
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotParameterizedCurve1.png\"/></p>
</html>"));
end parameterizedCurves;
