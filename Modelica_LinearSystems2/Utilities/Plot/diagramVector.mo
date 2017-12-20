within Modelica_LinearSystems2.Utilities.Plot;
function diagramVector "Plot several diagrams in vector layout"
  extends Modelica.Icons.Function;

  input Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram[:]
    "Properties of a set of diagrams (vector layout)"
    annotation(Dialog);
  input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
    Modelica_LinearSystems2.Utilities.Plot.Records.Device()
    "Properties of device where figure is shown" annotation(Dialog);

  import Modelica_LinearSystems2.Utilities.Plot.Internal;
protected
  Real mmToPixel= device.windowResolution/25.4;
  Integer position[4];
  Integer style;
  Integer nDiagrams=size(diagram,1);
  Integer figureHeight;
  Integer height[nDiagrams];
  Boolean OK;
  Integer id;

  function round "Round to nearest Integer"
    input Real r;
    output Integer i;
  algorithm
    i :=if r > 0 then integer(floor(r + 0.5)) else integer(ceil(r - 0.5));
  end round;

algorithm
  id := -1;

  // Determine total size of window
  figureHeight := 0;
  for i in 1:nDiagrams loop
    height[i]    := round(diagram[i].heightRatio*device.diagramWidth*mmToPixel);
    figureHeight := figureHeight + height[i];
  end for;

  position := {round(device.xTopLeft*mmToPixel),
               round(device.yTopLeft*mmToPixel),
               round(device.diagramWidth*mmToPixel),
               figureHeight};

  // Determine position and size of sub-figures
  for i in 1:nDiagrams loop
    if i > 1 then
      position[4] := height[i];
    end if;

    id:= createPlot(id=id,
                    position=position,
                    erase=true,
                    autoscale=true,
                    autoerase=false,
                    subPlot=i,
                    heading=diagram[i].heading,
                    grid=diagram[i].grid,
                    logX=diagram[i].logX,
                    logY=diagram[i].logY,
                    bottomTitle=diagram[i].xLabel,
                    leftTitle=diagram[i].yLabel,
                    color=device.autoLineColor,
                    legend=diagram[i].legend,
                    legendHorizontal=diagram[i].legendHorizontal,
                    legendFrame=diagram[i].legendFrame,
                    legendLocation=Integer(diagram[i].legendLocation));

    for j in 1:size(diagram[i].curve,1) loop
      if diagram[i].curve[j].autoLine then
        OK :=plotArray(diagram[i].curve[j].x,
                       diagram[i].curve[j].y,
                       legend=diagram[i].curve[j].legend,
                       id=id,
                       erase=false);
      else
        OK :=plotArray(diagram[i].curve[j].x,
                       diagram[i].curve[j].y,
                       legend=diagram[i].curve[j].legend,
                       color=diagram[i].curve[j].lineColor,
                       pattern=Internal.convertToDymolaPattern(diagram[i].curve[j].linePattern),
                       marker=Internal.convertToDymolaMarker(diagram[i].curve[j].lineSymbol),
                       thickness=diagram[i].curve[j].lineThickness,
                       id=id,
                       erase=false);
      end if;
    end for;
  end for;

  annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This function plots a set of 2-dimensional curves in a set of diagrams
using a vector layout. For an overview, see the documentation of package
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot\">Modelica_LinearSystems2.Utilities.Plot</a>.
</p>

<h4>Example</h4>
<p>
See an <a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot.Examples.showSinesInVectorDiagrams\">example</a>
for possible usage of this function:
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/showSinesInVectorDiagrams.png\">
</blockquote>
</html>"));
end diagramVector;
