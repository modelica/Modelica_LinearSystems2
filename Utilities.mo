within Modelica_LinearSystems2;
package Utilities "Functions that shall be included in Modelica.Utilities"
    extends Modelica.Icons.Package;
package Plot "Functions for generation of 2D-plots"
    extends Modelica.Icons.Package;
  package Examples "Demonstrate the usage of the plot functions"
     extends Modelica.Icons.ExamplesPackage;
    function plotSine "Plot a sine function in one diagram"
       input Modelica.SIunits.Frequency freqHz = 2 "Frequency of sine wave";
       input Modelica.SIunits.Damping damping = 0.8
          "Damping coefficient of sine wave";
      protected
       Integer nX = 500;
       Integer nPeriod = 5;
       Real x[nX];
       Real y[nX];
    algorithm
       (x,y) :=Utilities.dampedSine(freqHz, damping, nPeriod, nX);

      Modelica_LinearSystems2.Utilities.Plot.diagram(
        Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
            curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x,
              y=y,
              legend="torque1")},
            heading="Bearing friction torque",
            xLabel="w [rad/s]",
            yLabel="[N.m]"));
        annotation(interactive=true);
    end plotSine;

    function plotTwoSine "Plot two sine functions in one diagram"
       input Modelica.SIunits.Frequency freqHz1 = 2 "Frequency of sine wave 1";
       input Modelica.SIunits.Frequency freqHz2 = 3 "Frequency of sine wave 2";
       input Modelica.SIunits.Damping damping1 = 0.8
          "Damping coefficient of sine wave 12";
       input Modelica.SIunits.Damping damping2 = 0.1
          "Damping coefficient of sine wave 2";

      protected
       Integer nX = 500;
       Integer nSymbol = 20;
       Integer nPeriod = 5;
       Real x1[nX];
       Real y1[nX];
       Real x2[nX];
       Real y2[nX];
    algorithm
       (x1,y1) :=Utilities.dampedSine(freqHz1, damping1, nPeriod, nX);
       (x2,y2) :=Utilities.dampedSine(freqHz2, damping2, nPeriod, nX);

      Modelica_LinearSystems2.Utilities.Plot.diagram(
        Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
            curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x1,
              y=y1,
              legend="torque1"),
          Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x2,
              y=y2,
              legend="torque2")},
            heading="Bearing friction torques",
            xLabel="w [rad/s]",
            yLabel="[N.m]"));
        annotation(interactive=true);
    end plotTwoSine;

    function showSinesInVectorDiagrams
        "Plot sine functions in several diagrams in vector layout"
       input Modelica.SIunits.Frequency freqHz1 = 2 "Frequency of sine wave 1";
       input Modelica.SIunits.Frequency freqHz2 = 3 "Frequency of sine wave 2";
       input Modelica.SIunits.Damping damping1 = 0.8
          "Damping coefficient of sine wave 12";
       input Modelica.SIunits.Damping damping2 = 0.1
          "Damping coefficient of sine wave 2";

      protected
       Integer nX = 500;
       Integer nSymbol = 50;
       Integer nPeriod = 5;
       Real x1[nX];
       Real y1[nX];
       Real x2[nX];
       Real y2[nX];
       Real x3[nSymbol];
       Real y3[nSymbol];
    algorithm
       (x1,y1) :=Utilities.dampedSine(freqHz1, damping1, nPeriod, nX);
       (x2,y2) :=Utilities.dampedSine(freqHz2, damping2, nPeriod, nX);
       (x3,y3) :=Utilities.dampedSine(freqHz1, damping1, nPeriod, nSymbol);

      Modelica_LinearSystems2.Utilities.Plot.diagramVector({
        Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
            curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x1,
              y=y1,
              legend="torque1"),
          Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x2,
              y=y2,
              legend="torque2")},
            heading="Bearing friction torques",
            yLabel="[N.m]"),
        Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
            curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x1,
              y=1.1*y1,
              legend="torque3"),
          Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x2,
              y=1.2*y2,
              legend="torque4"),
          Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
              x=x3,
              y=y3,
              autoLine=false,
              linePattern=Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.None,
              lineSymbol=Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.Circle)},
            xLabel="w [rad/s]",
            yLabel="[N.m]")});
        annotation(interactive=true);
    end showSinesInVectorDiagrams;

    function showLegendStyles
        "Show several vector-diagram plots that demonstrate the various legend options"
      import Modelica_LinearSystems2.Utilities.Plot.Records;
      import Modelica_LinearSystems2.Utilities.Plot.Types;

      protected
       Real x[2]={0,1};
       Real y1[2]={0,0.5};
       Real y2[2]={0.2,0.7};
       Real y3[2]={0.4,0.9};
       Records.Curve curves[3]=
                           {Records.Curve(
                              x=x,
                              y=y1,
                              legend="curve 1"),
                            Records.Curve(
                              x=x,
                              y=y2,
                              legend="curve 2"),
                            Records.Curve(
                              x=x,
                              y=y3,
                              legend="curve 3")};
    algorithm
      Modelica_LinearSystems2.Utilities.Plot.diagramVector({Records.Diagram(
            curve=curves,
            heading="LegendLocation.Above",
            legendLocation=Types.LegendLocation.Above),Records.Diagram(
            curve=curves,
            heading="LegendLocation.Above, legendHorizontal=false",
            legendHorizontal=false,
            legendLocation=Types.LegendLocation.Above),Records.Diagram(
            curve=curves,
            heading="LegendLocation.Right",
            legendLocation=Types.LegendLocation.Right),Records.Diagram(
            curve=curves,
            heading="LegendLocation.Below",
            legendLocation=Types.LegendLocation.Below)}, Records.Device(
            xTopLeft=0,
            yTopLeft=0,
            diagramWidth=100));

      Modelica_LinearSystems2.Utilities.Plot.diagramVector({Records.Diagram(
            curve=curves,
            heading="LegendLocation.TopLeft",
            legendLocation=Types.LegendLocation.TopLeft),Records.Diagram(
            curve=curves,
            heading="LegendLocation.TopRight, legendHorizontal=false",
            legendHorizontal=false,
            legendLocation=Types.LegendLocation.TopRight),Records.Diagram(
            curve=curves,
            heading="LegendLocation.BottomLeft",
            legendLocation=Types.LegendLocation.BottomLeft),Records.Diagram(
            curve=curves,
            heading="LegendLocation.BottomRight",
            legendLocation=Types.LegendLocation.BottomRight)}, Records.Device(
            xTopLeft=105,
            yTopLeft=0,
            diagramWidth=100));

          annotation(interactive=true);
    end showLegendStyles;

    function showMatrixDiagrams
        "Demonstrate the layout of diagrams in matrix layout"
      import Modelica_LinearSystems2.Utilities.Plot.Records;
      import Modelica_LinearSystems2.Utilities.Plot.Types;

      protected
       Real x[2]={0,1};
       Real y1[2]={0,0.5};
       Real y2[2]={0.2,0.7};
       Real y3[2]={0.4,0.9};
       Records.Curve curves[3]=
                           {Records.Curve(
                              x=x,
                              y=y1,
                              legend="curve 1"),
                            Records.Curve(
                              x=x,
                              y=y2,
                              legend="curve 2"),
                            Records.Curve(
                              x=x,
                              y=y3,
                              legend="curve 3")};
    algorithm
      Modelica_LinearSystems2.Utilities.Plot.diagramMatrix(fill(
            Records.Diagram(curve=curves),
            5,
            3), Records.Device(diagramWidth=80));

      annotation(interactive=true);
    end showMatrixDiagrams;

    package Utilities
        "Utility functions (usually not of interest for the user)"
        extends Modelica.Icons.Package;
      function dampedSine "Return a damped sine chracteristic"
         input Modelica.SIunits.Frequency freqHz "Frequency of sine wave";
         input Modelica.SIunits.Damping damping
            "Damping coefficient of sine wave";
         input Integer nPeriod=5 "Number of periods to show";
         input Integer nPoints(min=2)=500 "Number of points";
         output Real x[nPoints];
         output Real y[nPoints];
        protected
       Real xEnd = 1/freqHz*nPeriod;
      algorithm
         x :=0:xEnd/(nPoints - 1):xEnd;
         for i in 1:size(x,1) loop
            y[i] :=Modelica.Math.exp(-x[i]*damping)*
                   Modelica.Math.sin(2*Modelica.Constants.pi*freqHz*x[i]);
         end for;
      end dampedSine;
    end Utilities;
  end Examples;

  function diagram "Plot one diagram"
    import Modelica_LinearSystems2.Utilities.Plot.Types;
    import Modelica.Utilities.Streams.*;
     input Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram
        "Diagram to be shown"
                            annotation(Dialog);
     input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
        Modelica_LinearSystems2.Utilities.Plot.Records.Device()
        "Properties of device where figure is shown"
                                                   annotation(Dialog);
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

     for i in 1:nCurves loop
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
      end for;

    annotation (interactive = true, Documentation(info="<html>
<p>
This function plots a set of 2-dimensional curves in a diagram.
For an overview, see the documentation of package
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot\">Modelica_LinearSystems2.Utilities.Plot</a>.
</p>
 
</html>"));
  end diagram;

  function diagramVector "Plot several diagrams in vector layout"
     input Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram[:]
        "Properties of a set of diagrams (vector layout)"                                                     annotation(Dialog);
     input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
        Modelica_LinearSystems2.Utilities.Plot.Records.Device()
        "Properties of device where figure is shown"
                                                   annotation(Dialog);

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
                     legendLocation=diagram[i].legendLocation);

         for j in 1:size(diagram[i].curve,1) loop
            if diagram[i].curve[j].autoLine or
               diagram[i].curve[j].lineSymbol==Types.PointSymbol.None then
               style :=0;
            elseif diagram[i].curve[j].linePattern==Types.LinePattern.None then
               style :=-(diagram[i].curve[j].lineSymbol - 1);
            else
               style :=diagram[i].curve[j].lineSymbol - 1;
            end if;

            OK :=plotArray(diagram[i].curve[j].x,
                           diagram[i].curve[j].y,
                           legend=diagram[i].curve[j].legend,
                           style=style,
                           id=id);
         end for;
      end for;

    annotation (interactive=true, Documentation(info="<html>
 
<p>
This function plots a set of 2-dimensional curves in a set of diagrams
using a vector layout. For an overview, see the documentation of package
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot\">Modelica_LinearSystems2.Utilities.Plot</a>.
</p>
 
</html>"));
  end diagramVector;

  function diagramMatrix "Plot several diagrams in matrix layout"
     input Modelica_LinearSystems2.Utilities.Plot.Records.Diagram diagram[:,:]
        "Properties of a set of diagrams (matrix layout)"                                                     annotation(Dialog);
     input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
        Modelica_LinearSystems2.Utilities.Plot.Records.Device()
        "Properties of device where figure is shown"
                                                   annotation(Dialog);
    protected
    Modelica_LinearSystems2.Utilities.Plot.Records.Device device2=device;

  algorithm
     for i in 1:size(diagram,2) loop
        device2.xTopLeft :=device.xTopLeft + (i - 1)*device.diagramWidth;
      Modelica_LinearSystems2.Utilities.Plot.diagramVector(diagram[:, i], device2);
     end for;

    annotation (interactive=true, Documentation(info="<html>
 
<p>
This function plots a set of 2-dimensional curves in a set of diagrams
using a matrix layout. For an overview, see the documentation of package
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Plot\">Modelica_LinearSystems2.Utilities.Plot</a>.
</p>
 
</html>"));
  end diagramMatrix;

  package Types "Types used for the plotting functions"
    extends Modelica.Icons.Package;
    type DrawingUnit_mm "Drawing unit specifing the nominal size in [mm]"
       extends Modelica.Icons.TypeReal(final quantity="Length", final unit="mm");
      annotation (Documentation(info="<html>
<p>
All size information for plotting, such as width or height
of a window or the thickness of a line, are defined by type <b>DrawingUnit_mm</b>. The 
recommended interpretation is that the DrawingUnit is the unscaled size in a document or
on printer in [mm]. For example, if the width of a diagram is 120,
and the diagram is pasted into a Word or PowerPoint document, 
then the width of the diagram in the document is 120 mm.
</p>
 
</html>"));
    end DrawingUnit_mm;

  type ImageResolution_dpi
        "Resolution of image in pixel per inch (screen) or dots per inch (printer)"
     extends Modelica.Icons.TypeReal(final quantity="ImageResolution");

    annotation (Documentation(info="<html>
<p>\"ImageResolution_dpi\" defines the mapping of a length coordinate to the resolution of the output device. The resolution [dpi] is defined as \"dots-per-inch\" and therefore a length L_mm defined in [mm] is mapped to a length L_dot in dots (or pixel) with the formula: </p>

<pre>   L_dot = round(ImageResolution_dpi/25.4 * L_mm)
</pre>

<p>where function round(..) rounds to the nearest integer. Typical values are \"96 dpi\" (for screen) or \"600 dpi\" for printer. For example if an \"ImageResolution = 96 dpi\" shall be used for a screen, then 1 mm is mapped to 4 pixel. </p>
</html>"));
  end ImageResolution_dpi;

  type LinePattern = enumeration(
          None,
          Solid) "Choices for line pattern";

  type PointSymbol = enumeration(
          None "no symbol",
          Cross "x",
          Circle "o",
          Square "square") "Choices for point symbol";

  type LegendLocation = enumeration(
          Above "Above diagram",
          Right "Right of diagram",
          Below "Below of diagram",
          TopLeft "Top left corner of diagram",
          TopRight "Top right corner of diagram",
          BottomLeft "Bottom left corner of diagram",
          BottomRight "Bottom right corner of diagram")
        "Choices for legend location";

  end Types;

  package Records "Records used to define the function interfaces"
   extends Modelica.Icons.Package;
    record Diagram
        "Properties of a diagram in a figure containing one or more curves"
      extends Modelica.Icons.Record;

      Modelica_LinearSystems2.Utilities.Plot.Records.Curve curve[:]
          "Properties of the curves in one diagram of the figure"  annotation(Dialog);

      String heading="" "Heading displayed above diagram" annotation(Dialog);
      Real heightRatio = 0.45 "Height of diagram = heightRatio*diagramWidth" annotation(Dialog);
      Boolean grid=true "= true, if grid is shown" annotation(Dialog,  choices(__Dymola_checkBox=true));

      /* group "Axes" (Axes properties) */
      String xLabel="" "String displayed at horizontal axis" annotation(Dialog(group="Axes"));
      String yLabel="" "String displayed at vertical axis" annotation(Dialog(group="Axes"));
      Boolean logX = false "= true, if logarithmic scale of x-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
      Boolean logY = false "= true, if logarithmic scale of y-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
      Boolean uniformScaling = false
          "= true, if same vertical and horizontal axis increment"
          annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));

      /* group "Legend" (Legend properties) */
      Boolean legend = true "= true, if legend is shown" annotation(Dialog(group="Legend"),choices(__Dymola_checkBox=true));
      Boolean legendFrame=false "= true, if frame around legend"
            annotation(Dialog(group="Legend"),   choices(__Dymola_checkBox=true));
      Boolean legendHorizontal=true
          "= true, if horizontal legend (provided it is meaningful)"
            annotation(Dialog(group="Legend"),choices(__Dymola_checkBox=true));
      Modelica_LinearSystems2.Utilities.Plot.Types.LegendLocation
          legendLocation=
                       Modelica_LinearSystems2.Utilities.Plot.Types.LegendLocation.Above
          "Legend placement"
                           annotation(Dialog(group="Legend"));

    end Diagram;

    record Curve "Properties of a curve (displayed in a diagram)"
      extends Modelica.Icons.Record;

       Real x[:] "x-values of curve" annotation(Dialog);
       Real y[:] "y-values of curve" annotation(Dialog);
       String legend="" "Legend text of curve" annotation(Dialog);

       Boolean autoLine = true "= true, if automatic line properties of curve"
         annotation(Dialog,  choices(__Dymola_checkBox=true));

       Integer lineColor[3]={0,0,255} "Color of curve as rgb values"
         annotation(Dialog(group="If autoLine = false (otherwise ignored)",__Dymola_colorSelector, __Dymola_treeView=false));

       Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern linePattern=
          Modelica_LinearSystems2.Utilities.Plot.Types.LinePattern.Solid
          "Line pattern of curve"
                                annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
       Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol lineSymbol=
          Modelica_LinearSystems2.Utilities.Plot.Types.PointSymbol.None
          "Symbol for points on curve"
                                     annotation(Dialog(group="If autoLine = false (otherwise ignored)"));

    /*
   Modelica_LinearSystems2.Utilities.Plot.Types.LineThickness_mm lineThickness=0.25 
    "Line thickness of curve" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
   Modelica_LinearSystems2.Utilities.Plot.Types.LineThickness_mm lineSymbolSize=3 
    "Symbol size" annotation(Dialog(group="If autoLine = false (otherwise ignored)"));
*/

    end Curve;

    record Device "Properties of a device"
      extends Modelica.Icons.Record;

      Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm xTopLeft=0
          "Horizontal position of top left figure corner if applicable (e.g. window)"
        annotation(Dialog);
      Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm yTopLeft=0
          "Vertical position of top left figure corner if applicable (e.g. window)"
        annotation(Dialog);
      Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm diagramWidth=140
          "Width of diagram"
                           annotation(Dialog);
      Modelica_LinearSystems2.Utilities.Plot.Types.ImageResolution_dpi
          windowResolution=
                         96
          "[dpi] Image resolution in window if applicable (e.g. unscaled window)"
                                                                                  annotation(Dialog);

      Boolean autoLineColor = true
          "if automatic line properties: distinguish curves by color otherwise by line style"
           annotation(Dialog,choices(__Dymola_checkBox=true));
    end Device;

    record DefaultDiagram
        "Diagram with no curves (might be used as default diagram)"
      extends Diagram(curve=fill(Curve(
                                 x=fill(0.0,0),y=fill(0.0,0)),0));
    end DefaultDiagram;
  end Records;
  annotation (Documentation(info="<html>
<p>
This package provides functions to plot curves in two dimensions. Here is a short overview:
</p>

<p>
A figure consists of a <u>set of diagrams</u>. Different functions are provided
to either plot one diagram (Plot.diagram) to plot several diagrams under each other
(Plot.vectorDiagrams) or to plot several diagrams in matrix layout
(Plot.matrixDiagrams).
</p>

<p>
Every diagram can have a set of <u>curves</u>. Every diagram has the same
width, defined by <u>diagramWidth</u>. The height of a diagram is defined
by variable <u>heightRatio</u> 
(diagram height in row j = diagram[j].heightRatio*diagramWidth). 
Several curves can be displayed in one diagram.
</p>
 
</html>"));
end Plot;

package Import "Functions to import data in a Modelica environment"
 extends Modelica.Icons.Package;
  package Examples "Demonstrate the usage of the Import functions"
      extends Modelica.Icons.ExamplesPackage;
    function linearizeDoublePendulum "Linearize double pendulum"
      output Real A[:,:] "A-matrix";
      output Real B[:,:] "B-matrix";
      output Real C[:,:] "C-matrix";
      output Real D[:,:] "D-matrix";
      output String inputNames[:] "Modelica names of inputs";
      output String outputNames[:] "Modelica names of outputs";
      output String stateNames[:] "Modelica names of states";
    algorithm
      (A,B,C,D,inputNames,outputNames,stateNames) :=
        Modelica_LinearSystems2.Utilities.Import.linearize("Modelica_LinearSystems2.Utilities.Import.Examples.Utilities.DoublePendulum",
        1.0);
      annotation(interactive=true);
    end linearizeDoublePendulum;

    package Utilities
        extends Modelica.Icons.Library;
      model DoublePendulum "double pendulum system"

        parameter Modelica.SIunits.Mass m_trolley = 5;
        parameter Modelica.SIunits.Mass m_load = 20;
        parameter Modelica.SIunits.Length length = 2;
        parameter Modelica.SIunits.Angle phi1_start = -80.0/180*pi;
        parameter Modelica.SIunits.Angle phi2_start = 10;
        parameter Modelica.SIunits.AngularVelocity w1_start = 0.0;
        parameter Modelica.SIunits.AngularVelocity w2_start = 0.0;

        constant Real pi = Modelica.Constants.pi;

        inner Modelica.Mechanics.MultiBody.World world(gravityType=Modelica.Mechanics.MultiBody.Types.GravityTypes.
              UniformGravity, animateWorld=false)
                              annotation (Placement(transformation(extent={{-140,-80},
                  {-120,-60}},
                            rotation=0)));
        Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
          annotation (Placement(transformation(extent={{-96,0},{-76,20}})));
        Modelica.Mechanics.Translational.Components.Damper damper1(d=0)
          annotation (Placement(transformation(extent={{-96,14},{-76,34}})));
        Modelica.Mechanics.MultiBody.Joints.Revolute rev(n={0,0,1},useAxisFlange=true,
          phi(fixed=true, start=phi1_start),
          w(fixed=true, start=w1_start))
                                     annotation (Placement(transformation(extent={{-30,0},
                  {-10,20}},      rotation=0)));
        Modelica.Mechanics.Rotational.Components.Damper damper(d=0)
          annotation (Placement(transformation(extent={{-22,40},{-2,60}},rotation=0)));
        Modelica.Mechanics.MultiBody.Parts.Body body(
          m=m_load,
          r_CM={0,0,0},
          specularCoefficient=4*world.defaultSpecularCoefficient,
          sphereDiameter=1.5*world.defaultBodyDiameter)
          annotation (Placement(transformation(extent={{78,0},{98,20}}, rotation=0)));
        Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
          shapeType="box",
          animateSphere=true,
          m=m_trolley,
          sphereDiameter=world.defaultBodyDiameter)
          annotation (Placement(transformation(extent={{-58,0},{-38,20}})));
        Modelica.Mechanics.Translational.Sources.Force force
          annotation (Placement(transformation(extent={{-98,34},{-78,54}})));
        Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles
          annotation (Placement(transformation(extent={{-30,-30},{-10,-10}})));
        Modelica.Mechanics.MultiBody.Sensors.RelativeVelocity relativeVelocity
          annotation (Placement(transformation(extent={{-96,-30},{-76,-10}})));
        Modelica.Mechanics.MultiBody.Sensors.RelativePosition relativePosition
          annotation (Placement(transformation(extent={{-96,-60},{-76,-40}})));
        Modelica.Blocks.Interfaces.RealInput u
          annotation (Placement(transformation(extent={{-190,-20},{-150,20}})));
        Modelica.Blocks.Interfaces.RealOutput s
          annotation (Placement(transformation(extent={{150,90},{170,110}})));
        Modelica.Blocks.Interfaces.RealOutput v
          annotation (Placement(transformation(extent={{150,50},{170,70}})));
       Modelica.Blocks.Interfaces.RealOutput phi
          annotation (Placement(transformation(extent={{150,10},{170,30}})));
        Modelica.Blocks.Interfaces.RealOutput w
          annotation (Placement(transformation(extent={{150,-30},{170,-10}})));
        Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
            relativeAngularVelocity
          annotation (Placement(transformation(extent={{-30,-60},{-10,-40}})));

        Modelica.Blocks.Sources.Constant const(k=0.5*Modelica.Constants.pi)
          annotation (Placement(transformation(extent={{94,-22},{106,-10}})));
        Modelica.Blocks.Math.Add add
          annotation (Placement(transformation(extent={{116,-10},{136,10}})));
        Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(
          phi(fixed=true, start=phi2_start),
          w(fixed=true, start=w2_start),
          cylinderDiameter=3*world.defaultJointWidth,
          cylinderColor={0,0,200})                             annotation (Placement(transformation(extent={{24,0},{
                  44,20}}, rotation=0)));
        Modelica.Mechanics.MultiBody.Sensors.RelativeAngles relativeAngles1
          annotation (Placement(transformation(extent={{24,-30},{44,-10}})));
        Modelica.Mechanics.MultiBody.Sensors.RelativeAngularVelocity
            relativeAngularVelocity1
          annotation (Placement(transformation(extent={{24,-60},{44,-40}})));
       Modelica.Blocks.Interfaces.RealOutput phi1
          annotation (Placement(transformation(extent={{150,-70},{170,-50}})));
        Modelica.Blocks.Interfaces.RealOutput w1
          annotation (Placement(transformation(extent={{150,-110},{170,-90}})));
        Modelica.Blocks.Math.Add add1
          annotation (Placement(transformation(extent={{88,-50},{108,-30}})));
        Modelica.Blocks.Sources.Constant const1(k=0)
          annotation (Placement(transformation(extent={{66,-62},{78,-50}})));
        Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
          r={length/2,0,0},
          specularCoefficient=0.7,
          color={0,0,0},
          diameter=0.05,
          density=900)
          annotation (Placement(transformation(extent={{-4,0},{16,20}})));
        Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
          r={length/2,0,0},
          specularCoefficient=0.7,
          color={0,0,0},
          diameter=0.05,
          density=900)
          annotation (Placement(transformation(extent={{52,0},{72,20}})));
      equation
        connect(damper.flange_b, rev.axis) annotation (Line(points={{-2,50},{0,50},{0,
                24},{0,20},{-20,20}},   color={0,0,0}));
        connect(rev.support, damper.flange_a) annotation (Line(points={{-26,20},{-26,
                26},{-36,26},{-36,50},{-22,50}}, color={0,0,0}));
        connect(bodyShape.frame_b, rev.frame_a) annotation (Line(
            points={{-38,10},{-30,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(prismatic.frame_a, world.frame_b) annotation (Line(
            points={{-96,10},{-110,10},{-110,-70},{-120,-70}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(force.flange, prismatic.axis) annotation (Line(
            points={{-78,44},{-78,16}},
            color={0,127,0},
            smooth=Smooth.None));
        connect(damper1.flange_a, prismatic.support) annotation (Line(
            points={{-96,24},{-96,16},{-90,16}},
            color={0,127,0},
            smooth=Smooth.None));
        connect(damper1.flange_b, prismatic.axis) annotation (Line(
            points={{-76,24},{-78,24},{-78,16}},
            color={0,127,0},
            smooth=Smooth.None));
        connect(prismatic.frame_b, bodyShape.frame_a) annotation (Line(
            points={{-76,10},{-58,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeVelocity.frame_b, prismatic.frame_b) annotation (Line(
            points={{-76,-20},{-76,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeVelocity.frame_a, prismatic.frame_a) annotation (Line(
            points={{-96,-20},{-96,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativePosition.frame_b, relativeVelocity.frame_b) annotation (Line(
            points={{-76,-50},{-76,-20}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativePosition.frame_a, relativeVelocity.frame_a) annotation (Line(
            points={{-96,-50},{-96,-20}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngles.frame_b, rev.frame_b) annotation (Line(
            points={{-10,-20},{-10,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngles.frame_a, rev.frame_a) annotation (Line(
            points={{-30,-20},{-30,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(u, force.f) annotation (Line(
            points={{-170,0},{-136,0},{-136,44},{-100,44}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(relativeAngularVelocity.frame_a, relativeAngles.frame_a) annotation (
            Line(
            points={{-30,-50},{-30,-20}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngularVelocity.frame_b, relativeAngles.frame_b) annotation (
            Line(
            points={{-10,-50},{-10,-20}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngularVelocity.w_rel[3], w) annotation (Line(
            points={{-20,-60.3333},{-20,-66},{120,-66},{120,-20},{160,-20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(relativeVelocity.v_rel[1], v) annotation (Line(
            points={{-86,-31.6667},{-104,-31.6667},{-104,-32},{-118,-32},{-118,
                  62},{42,62},{42,60},{160,60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(relativePosition.r_rel[1], s) annotation (Line(
            points={{-86,-61.6667},{-104,-61.6667},{-104,-58},{-120,-58},{-120,
                  100},{160,100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(phi, phi) annotation (Line(
            points={{160,20},{160,20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add.y, phi) annotation (Line(
            points={{137,0},{148,0},{148,20},{160,20}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(const.y, add.u2) annotation (Line(
            points={{106.6,-16},{110,-16},{110,-6},{114,-6}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add.u1, relativeAngles.angles[3]) annotation (Line(
            points={{114,6},{108,6},{108,-4},{58,-4},{58,-36},{-20,-36},{-20,
                  -30.3333}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(relativeAngles1.frame_a, revolute2.frame_a) annotation (Line(
            points={{24,-20},{24,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngles1.frame_b, revolute2.frame_b) annotation (Line(
            points={{44,-20},{44,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngles1.frame_a, relativeAngularVelocity1.frame_a)
          annotation (Line(
            points={{24,-20},{24,-50}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(relativeAngularVelocity1.frame_b, relativeAngles1.frame_b)
          annotation (Line(
            points={{44,-50},{44,-20}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(const1.y, add1.u2)
                                 annotation (Line(
            points={{78.6,-56},{82,-56},{82,-46},{86,-46}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add1.u1, relativeAngles1.angles[3]) annotation (Line(
            points={{86,-34},{60,-34},{60,-30.3333},{34,-30.3333}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(add1.y, phi1) annotation (Line(
            points={{109,-40},{136,-40},{136,-60},{160,-60}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(relativeAngularVelocity1.w_rel[3], w1) annotation (Line(
            points={{34,-60.3333},{36,-60.3333},{36,-100},{160,-100}},
            color={0,0,127},
            smooth=Smooth.None));
        connect(bodyCylinder1.frame_b, body.frame_a) annotation (Line(
            points={{72,10},{78,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (Line(
            points={{52,10},{44,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(bodyCylinder.frame_b, revolute2.frame_a) annotation (Line(
            points={{16,10},{24,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        connect(bodyCylinder.frame_a, rev.frame_b) annotation (Line(
            points={{-4,10},{-10,10}},
            color={95,95,95},
            thickness=0.5,
            smooth=Smooth.None));
        annotation (
          experiment(
            StartTime=1,
            StopTime=10,
            Algorithm="Dassl"),
          Diagram(coordinateSystem(
              preserveAspectRatio=true,
              extent={{-150,-100},{150,100}},
              grid={2,2}), graphics),
          Documentation(info="<html>
 
<p>
Model of a simple double pendulum system. 
It is used to demonstrate linearization
</p>
 
</html>"),experimentSetupOutput,
          Icon(coordinateSystem(preserveAspectRatio=true, extent={{-150,-100},{
                    150,100}}), graphics={
                Rectangle(
                  extent={{-150,122},{150,-120}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid),
                Rectangle(
                  extent={{-82,22},{82,18}},
                  lineColor={0,0,255},
                  fillPattern=FillPattern.Forward),
                Rectangle(extent={{-44,54},{0,28}}, lineColor={0,0,0}),
                Ellipse(
                  extent={{-40,34},{-28,22}},
                  lineColor={0,0,0},
                  fillPattern=FillPattern.Solid,
                  fillColor={255,255,255},
                  lineThickness=0.5),
                Ellipse(
                  extent={{-16,34},{-4,22}},
                  lineColor={0,0,0},
                  fillColor={255,255,255},
                  fillPattern=FillPattern.Solid,
                  lineThickness=0.5),
                Line(
                  points={{-18,-16},{10,-62}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Ellipse(
                  extent={{4,-56},{20,-72}},
                  lineColor={0,0,0},
                  fillColor={0,255,255},
                  fillPattern=FillPattern.Solid),
                Ellipse(
                  extent={{-25,44},{-19,38}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid),
                Line(
                  points={{28,46},{4,46}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{34,40},{10,40}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Line(
                  points={{-22,40},{-18,-16}},
                  color={0,0,0},
                  smooth=Smooth.None),
                Ellipse(
                  extent={{-20,-15},{-14,-21}},
                  lineColor={0,0,0},
                  fillColor={95,95,95},
                  fillPattern=FillPattern.Solid)}));
      end DoublePendulum;
    end Utilities;
  end Examples;

  function linearize "Linearize a model after simulation up to a given time"
    input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
    input Modelica.SIunits.Time t_linearize= 0
        "Simulate until T_linearize and then linearize"
                                                      annotation(Dialog);

    protected
    String fileName="dslin";
    String fileName2=fileName+".mat";

    // Simulate until t_linearize and then linearize at this time instant
    Boolean OK1 = simulateModel(problem=modelName, startTime=0, stopTime=t_linearize);
    Boolean OK2 = importInitial("dsfinal.txt");
    Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=t_linearize, stopTime=t_linearize);

    // Read linear system from file
    Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
    Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
    Integer nx=integer(nxMat[1, 1]);
    Integer nu=ABCDsizes[2] - nx;
    Integer ny=ABCDsizes[1] - nx;
    Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
    String xuyName[nx + nu + ny]=readStringMatrix(fileName2, "xuyName", nx + nu + ny);
    public
    output Real A[nx,nx] =  ABCD[1:nx, 1:nx] "A-matrix";
    output Real B[nx,nu] =  ABCD[1:nx, nx + 1:nx + nu] "B-matrix";
    output Real C[ny,nx] =  ABCD[nx + 1:nx + ny, 1:nx] "C-matrix";
    output Real D[ny,nu] =  ABCD[nx + 1:nx + ny, nx + 1:nx + nu] "D-matrix";
    output String inputNames[nu] =  xuyName[nx + 1:nx + nu]
        "Modelica names of inputs";
    output String outputNames[ny] =  xuyName[nx + nu + 1:nx + nu + ny]
        "Modelica names of outputs";
    output String stateNames[nx] =  xuyName[1:nx] "Modelica names of states";
  algorithm

     annotation (interactive=true, Documentation(info="<html>
<p>This function initializes a Modelica model and then simulates the model with its default experiment options until time instant \"t_linearize\". If t_linearize=0, no simulation takes place (only initialization). At the simulation stop time, the model is linearized in such a form that </p>
<p><ul>
<li>all top-level signals with prefix \"input\" are treated as inputs <b>u</b>(t) of the model ,</li>
<li>all top-level signals with prefix \"output\" are treated as outputs <b>y</b>(t) of the model,</li>
<li>all variables that appear differentiated and that are selected as states at this time instant are treated as states <b>x</b> of the model.</li>
</ul></p>
<p>Formally, the non-linear hybrid differential-algebraic equation system is therefore treated as the following ordinary equation system at time instant t_linearize: </p>
<pre>    der(<b>x</b>) = <b>f</b>(<b>x</b>,<b>u</b>)</pre>
<pre>         <b>y</b> = <b>g</b>(<b>x</b>,<b>u</b>) </pre>
<p>Taylor series expansion (linearization) of this model around the simulation stop time t_linearize: </p>
<pre>   <b>u</b>0 = <b>u</b>(t_linearize)</pre>
<pre>   <b>y</b>0 = <b>y</b>(t_linearize)</pre>
<pre>   <b>x</b>0 = <b>x</b>(t_linearize) </pre>
<p>and neglecting higher order terms results in the following system: </p>
<pre>   der(<b>x</b>0+d<b>x</b>) = <b>f</b>(<b>x</b>0,<b>u</b>0) + der(<b>f</b>,<b>x</b>)*d<b>x</b> + der(<b>f</b>,<b>u</b>)*d<b>u</b></pre>
<pre>      <b>y</b>0 + d<b>y</b> = <b>g</b>(<b>x</b>0,<b>u</b>0) + der(<b>g</b>,<b>x</b>)*d<b>x</b> + der(<b>g</b>,<b>u</b>)*d<b>u</b></pre>
<p>where der(<b>f</b>,<b>x</b>) is the partial derivative of <b>f</b> with respect to <b>x</b>, and the partial derivatives are computed at the linearization point t_linearize. Re-ordering of terms gives (note <b>der</b>(<b>x</b>0) = <b>0</b>): </p>
<pre>   der(d<b>x</b>) = der(<b>f</b>,<b>x</b>)*d<b>x</b> + der(<b>f</b>,<b>u</b>)*d<b>u</b> + <b>f</b>(<b>x</b>0,<b>u</b>0)</pre>
<pre>        d<b>y</b> = der(<b>g</b>,<b>x</b>)*d<b>x</b> + der(<b>g</b>,<b>u</b>)*d<b>u</b> + (<b>g</b>(<b>x</b>0,<b>u</b>0) - <b>y</b>0)</pre>
<p>or </p>
<pre>   der(d<b>x</b>) = <b>A</b>*d<b>x</b> + <b>B</b>*d<b>u</b> + <b>f</b>0</pre>
<pre>        d<b>y</b> = <b>C</b>*d<b>x</b> + <b>D</b>*d<b>u</b></pre>
<p>This function returns the matrices <b>A</b>, <b>B</b>, <b>C</b>, <b>D</b> and assumes that the linearization point is a steady-state point of the simulation (i.e., <b>f</b>(<b>x</b>0,<b>u</b>0) = 0). Additionally, the full Modelica names of all inputs, outputs and states shall be returned if possible (default is to return empty name strings).</p>
</html>"));
  end linearize;
end Import;

end Utilities;
