within Modelica_LinearSystems2.Utilities.Plot;
function parameterizedCurves
  "Plot parametrized curve with one or more branches"
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.WorkInProgress.RootLocusOld.Types.MarkerStyles;

  input Modelica_LinearSystems2.Utilities.Plot.Records.ParametrizedCurves diagram
    "Parametrized curve data points" annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device
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

  id:= createPlot(id=-1,
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
                  legendLocation=diagram.legendLocation);

  // Plot parameterized curve
  if nProperties == 0 then
     plotParametricCurves(
                   x=diagram.X,
                   y=diagram.Y,
                   s=diagram.s,
                   xName=diagram.xName,
                   yName=diagram.yName,
                   sName=diagram.sName,
                   legends=diagram.legends,
                   id=  id,
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
     plotParametricCurves(
                   x=diagram.X,
                   y=diagram.Y,
                   s=diagram.s,
                   xName=diagram.xName,
                   yName=diagram.yName,
                   sName=diagram.sName,
                   legends=diagram.legends,
                   id=  id,
                   labelWithS=diagram.labelWithS,
                   colors=colors,
                   patterns=Internal.convertToDymolaPattern(patterns),
                   markers=Internal.convertToDymolaMarker(symbols),
                   thicknesses=thicknesses);
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

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
Plot.<b>rootLocus</b>(modelName, t_linearize, modelParam, simulationSetup, diagram)
</pre></blockquote>

<h4>Description</h4>
<p>
This function examines a root locus analysis of a selected Modelica model over 
variation of a certain system parameter. 
Note, only first parameter <code>modelParam[1]</code> is considered for the analysis. 
The parameter is varied equidistantly from minimum to maximum value. 
</p>

<h4>Example</h4>
<p>
Calling the function
</p>
<blockquote><pre>
Utilities.Plot.<b>rootLocus</b>(
  modelName = \"Modelica.Mechanics.Rotational.Examples.First\",
  t_linearize = 0,
  modelParam={
    Modelica_LinearSystems2.Records.ParameterVariation(
      Name=\"Jload\",
      Min=1,
      Max=6,
      nVar=10,
      Unit=\"kg.m2\")});
</pre></blockquote>
<p>
yields following diagram
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/></p>
</html>"));
end parameterizedCurves;
