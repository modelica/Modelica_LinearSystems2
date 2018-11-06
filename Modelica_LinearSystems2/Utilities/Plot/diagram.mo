within Modelica_LinearSystems2.Utilities.Plot;
function diagram "Plot one diagram"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Utilities.Plot.Types;
  import Modelica_LinearSystems2.Utilities.Plot.Internal;

  input Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram
    "Diagram to be shown" annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
    Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);
protected
  Real mmToPixel= device.windowResolution/25.4;
  Real position[4];
  Integer style;
  Integer nCurves=size(diagram.curve,1);
  Boolean OK;
  Integer id;
algorithm
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
                  thicknesses={0.25},
                  axes={1},
                  displayUnits={""});

  for i in 1:nCurves loop
    /*
    if diagram.curve[i].autoLine or
       diagram.curve[i].lineSymbol==Types.PointSymbol.None then
       style :=0;
    elseif diagram.curve[i].linePattern==Types.LinePattern.None then
       style :=-(diagram.curve[i].lineSymbol - 1);
    else
       style :=diagram.curve[i].lineSymbol - 1;
    end if;

    OK :=plotArray(diagram.curve[i].x,
                   diagram.curve[i].y,
                   legend=diagram.curve[i].legend,
                   style=style,
                   id=id);
    */

    if diagram.curve[i].autoLine then
      OK := DymolaCommands.Plot.plotArray(
        diagram.curve[i].x,
        diagram.curve[i].y,
        legend=diagram.curve[i].legend,
        id=id,
        erase=false);
    else
      OK := DymolaCommands.Plot.plotArray(
        diagram.curve[i].x,
        diagram.curve[i].y,
        legend=diagram.curve[i].legend,
        color=diagram.curve[i].lineColor,
        pattern=Internal.convertToDymolaPattern(diagram.curve[i].linePattern),
        marker=Internal.convertToDymolaMarker(diagram.curve[i].lineSymbol),
        thickness=diagram.curve[i].lineThickness,
        id=id,
        erase=false);
    end if;
  end for;

  annotation (
    Documentation(info="<html>
<p>
This function plots a set of 2-dimensional curves in a diagram.
For an overview, see the documentation of package
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot\">Modelica_LinearSystems2.Utilities.Plot</a>.
</p>

<h4>Example</h4>
<p>
See an <a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.plotTwoSineDifferentStyles\">example</a>
for possible usage of this function:
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/plotTwoSineDifferentStyles.png\">
</blockquote>
</html>"));
end diagram;
