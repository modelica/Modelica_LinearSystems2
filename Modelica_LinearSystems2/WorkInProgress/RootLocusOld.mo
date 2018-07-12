within Modelica_LinearSystems2.WorkInProgress;
package RootLocusOld
  function rootLocus
    "Plot root locus of nonlinear Modelica model by linearizing the model for variations of one model parameter"
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import
      Modelica_LinearSystems2.WorkInProgress.RootLocusOld.Types.MarkerStyles;

    input String modelName "Name of the Modelica model"
      annotation(Dialog(__Dymola_translatedModel(caption="Model to be linearized for the root locus")));

    input Modelica.SIunits.Time t_linearize= 0
      "Simulate until t_linearize and then linearize" annotation(Dialog);

  /*
  input Modelica.Mechanics.MultiBody.Interfaces.partialColorMap colorMap=
      Modelica.Mechanics.MultiBody.Visualizers.Colors.ColorMaps.jet
    "Color map function for parameter variations" annotation(choicesAllMatching=true);
*/

    input
      Modelica_LinearSystems2.WorkInProgress.RootLocusOld.ParameterVariation
                                                             modelParam[:]
      "Model parameter to be varied";

    input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
        Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
      "Simulation options it t_linearize > 0";

    input
      Modelica_LinearSystems2.WorkInProgress.RootLocusOld.RootLocusDiagramOld     diagram
      "Diagram properties of the root locus";

  protected
    Boolean ok "True, if call is ok";
    Real dp "Step of parameter equidistand grid";
    Real parValue "Value of parameter in a loop";
    Integer nVar = max(modelParam[1].nVar,2)
      "Number of parameter variations (at least 2 are necessary)";
    Real color[nVar,3] = Modelica.Mechanics.MultiBody.Visualizers.Colors.ColorMaps.jet(nVar)
      "Color of markers in a loop";
    Real mmToPixel= diagram.windowResolution/25.4;
    Real position[4];
    Real A[:,:] "System matrix of linearization";
    Real eigenValues[:,2] "Eigen values of A";
    Integer id;
    String heading;
  algorithm
    assert(size(modelParam,1) == 1, "Input vector modelParam must have length one, but has length "+ String(size(modelParam,1)));

    // Define heading
    if diagram.heading=="" then
       heading :="Root locus of " + modelName + " over " + modelParam[1].Name;
    else
       heading :=diagram.heading;
    end if;

    // Create diagram
    position := {diagram.xTopLeft,
                 diagram.yTopLeft,
                 diagram.diagramWidth,
                 diagram.heightRatio*diagram.diagramWidth}*mmToPixel;
    id:= createPlot(id=-1,
                    position=integer(position),
                    erase=true,
                    autoscale=true,
                    autoerase=false,
                    subPlot=1,
                    heading=heading,
                    grid=diagram.grid,
                    logX=diagram.logX,
                    logY=diagram.logY,
                    bottomTitle=diagram.xLabel,
                    leftTitle=diagram.yLabel,
                    legend=true,
                    legendHorizontal=false,
                    legendLocation=Modelica_LinearSystems2.Utilities.Plot.Types.LegendLocation.Right);

    // Initialize computation
    dp := (modelParam[1].Max-modelParam[1].Min)/(nVar-1);
    ok:=translateModel(modelName);
    Modelica.Utilities.Files.removeFile("dsin.mat");

    // Loop over the selected parameter range
    for i in 1:nVar loop
      // Linearize
      parValue := modelParam[1].Min+(i-1)*dp;
      ok :=SetVariable(modelParam[1].Name, parValue);
      A :=Modelica_LinearSystems2.WorkInProgress.RootLocusOld.linearize2(
              modelName,
              t_linearize=t_linearize,
              simulationSetup=simulationSetup);

      // Compute eigen values
      eigenValues :=Modelica.Math.Matrices.eigenValues(A);

      // Plot root loci
      ok :=plotArray(eigenValues[:,1],
                     eigenValues[:,2],
                     legend=if i==1 or i==nVar or mod(i,nVar/10) == 0 then String(parValue,significantDigits=3) else "",
                     color=integer(color[i, :]),
                     pattern=LinePattern.None,
                     marker=MarkerStyle.Square,
                     id=id,
                     erase=false);
    end for;
    Modelica.Utilities.Files.removeFile("dslin.mat");

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
  modelName = &quot;Modelica.Mechanics.Rotational.Examples.First&quot;,
  t_linearize = 0,
  modelParam={
    Modelica_LinearSystems2.Records.ParameterVariation(
      Name=&quot;Jload&quot;,
      Min=1,
      Max=6,
      nVar=10,
      Unit=&quot;kg.m2&quot;)});
</pre></blockquote>
<p>
yields following diagram
</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/></p>
</html>"));
  end rootLocus;

  record ParameterVariation
    "Define variation of one parameter in a given range and optionally select the parameter from a translated model"
    String Name "Name of parameter" annotation (Dialog);
    Real Min "Minimum value of parameter" annotation (Dialog);
    Real Max "Maximum value of parameter" annotation (Dialog);
    Integer nVar(min=2) = 10 "Number of parameter variations (min=2)" annotation (Dialog);
    String Unit="" "Unit of parameter" annotation (Dialog);
    annotation (
    Dialog(__Dymola_importDsin(onlyStart=true,
      fields(Name=initialName,
             Min=initialValue.minimum,
             Max=initialValue.maximum,
             Unit=initialValue.unit))),
    Icon(graphics={
          Rectangle(
            extent={{-100,-30},{100,-90}},
            lineColor={0,0,0},
            fillColor={175,175,175},
            fillPattern=FillPattern.Solid),
          Line(
            points={{-100,-30},{-100,44}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{100,-32},{100,44}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-100,40},{100,40}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{-60,68},{-100,40},{-60,12}},
            color={0,0,0},
            smooth=Smooth.None),
          Line(
            points={{60,68},{100,40},{60,12}},
            color={0,0,0},
            smooth=Smooth.None)}));
  end ParameterVariation;

  function rootLocusOfDriveOld
    "Plot the root locus of a drive with varying load"
  algorithm
    Modelica_LinearSystems2.WorkInProgress.RootLocusOld.rootLocus(
       "Modelica.Mechanics.Rotational.Examples.First",
       modelParam={Modelica_LinearSystems2.WorkInProgress.RootLocusOld.ParameterVariation(
                                                                      Name="Jload", Min=1, Max=20, nVar=30, Unit="kg.m2")});
      annotation(__Dymola_interactive=true, Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica.Mechanics.Rotational.Examples.First\">Rotational.Examples.First</a>
over the load inertia <b>Jload</b>:
</p>

<blockquote>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/rootLocusOfDrive.png\">
</blockquote>
</html>"));
  end rootLocusOfDriveOld;

  record RootLocusDiagramOld "Properties of a root locus diagram"
    extends Modelica.Icons.Record;

    String heading=""
      "Heading displayed above diagram (if empty, default heading)"                 annotation(Dialog);
    Real heightRatio = 0.8 "Height of diagram = heightRatio*diagramWidth" annotation(Dialog);
    Boolean grid=true "True, if grid is shown" annotation(Dialog,  choices(checkBox=true));

    /* group "Axes" (Axes properties) */
    String xLabel="Real part of eigenvalues"
      "String displayed at horizontal axis"                                        annotation(Dialog(group="Axes"));
    String yLabel="Imaginary part of eigenvalues"
      "String displayed at vertical axis"                                             annotation(Dialog(group="Axes"));
    Boolean logX = false "True, if logarithmic scale of x-axis" annotation(Dialog(group="Axes"),choices(checkBox=true));
    Boolean logY = false "True, if logarithmic scale of y-axis" annotation(Dialog(group="Axes"),choices(checkBox=true));
    Boolean uniformScaling = false
      "True, if same vertical and horizontal axis increment"
        annotation(Dialog(group="Axes"),choices(checkBox=true));

    Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm xTopLeft=0
      "Horizontal position of top left figure corner if applicable (e.g. window)"
      annotation(Dialog);
    Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm yTopLeft=0
      "Vertical position of top left figure corner if applicable (e.g. window)"
      annotation(Dialog);
    Modelica_LinearSystems2.Utilities.Plot.Types.DrawingUnit_mm diagramWidth=140
      "Width of diagram" annotation(Dialog);
    Modelica_LinearSystems2.Utilities.Plot.Types.ImageResolution_dpi windowResolution=96
      "[dpi] Image resolution in window if applicable (e.g. unscaled window)"   annotation(Dialog);

  end RootLocusDiagramOld;

  function linearize2
    "Linearize a model after simulation up to a given time and return only the A matrix"
    import Modelica.Utilities.Streams;

    input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
    input Modelica.SIunits.Time t_linearize= 0
      "Simulate until t_linearize and then linearize" annotation(Dialog);
    input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
        Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
      "Simulation options it t_linearize > 0" annotation(Dialog);

  protected
    String fileName="dslin";
    String fileName2=fileName+".mat";

    // Simulate until t_linearize and then linearize at this time instant
    Boolean OK1 = if t_linearize <= 0.0 then true else
                     simulateModel(problem=modelName, startTime=0, stopTime=t_linearize,
                                   method=simulationSetup.method,
                                   tolerance=simulationSetup.tolerance,
                                   fixedstepsize=simulationSetup.fixedStepSize);
    Boolean OK2 = if t_linearize <= 0.0 then true else importInitial("dsfinal.txt");
    Boolean OK3 = linearizeModel(problem=modelName, resultFile=fileName, startTime=t_linearize, stopTime=t_linearize);

    // Read linear system from file
    Integer ABCDsizes[2]=Streams.readMatrixSize(fileName2, "ABCD");
    Integer nx = integer(scalar(Streams.readRealMatrix(fileName2, "nx", 1, 1)));
    Integer nu=ABCDsizes[2] - nx;
    Integer ny=ABCDsizes[1] - nx;
    Real ABCD[nx + ny,nx + nu]=Streams.readRealMatrix(fileName2, "ABCD", nx + ny, nx + nu);
  public
    output Real A[nx,nx] =  ABCD[1:nx, 1:nx] "A-matrix";
  algorithm

     annotation (__Dymola_interactive=true, Documentation(info="<html>
<p>
This function is identical to function
<a href=\"modelica://Modelica_LinearSystems2.Utilities.Import.linearize\">linearize</a>
but returns only the A-matrix.
</p>
</html>"));
  end linearize2;

  package Types "Package of type definitions"
    extends Modelica.Icons.TypesPackage;

    type MarkerStyles = enumeration(
        Cross "Cross",
        Circle "Circle",
        FilledCircle "FilledCircle",
        Square "Square",
        FilledSquare "FilledSquare",
        TriangleDown "TriangleDown",
        TriangleUp "TriangleUp",
        Diamond "Diamond") "Style of marker of plotted curve";
    annotation (
      Documentation(info="<html>
<p>
This package contains type definitions used in the library. Generally,
the enumeration type is used to assign a unique choice of parameter
within a model.
</p>
</html>"));
  end Types;

  function plotEigenvalues
    "Calculate eigenvalues of matrix A and plot root locus"
    import
      Modelica_LinearSystems2.WorkInProgress.RootLocusOld.Types.MarkerStyles;

    input Real A[:,size(A, 1)] = [2,1,1;1,1,1;1,2,2] "Square matrix";
    input Boolean removePrevious=true
      "True, if all previous plots should be deleted"
      annotation (Dialog(group="Plot settings"));
    input Integer position[4]={5, 5, 600, 450} "Window Position"
      annotation (Dialog(group="Plot settings"));
    input String heading="Root locii" "Heading of plot"
      annotation (Dialog(group="Plot settings"));
    input Boolean useLegend = true "Use legend"
      annotation (Dialog(group="Plot settings"));
    input String legend="Eigenvalues" "Legend of plot"
      annotation (Dialog(group="Plot settings", enable=useLegend));
    input Boolean grid = false "Add grid"
      annotation (Dialog(group="Plot settings"));
    input MarkerStyles markerStyle=MarkerStyles.Cross "Style of marker"
      annotation (Dialog(group="Plot settings"));
    input Integer markerColor[3]={0,0,255} "Color of marker"
      annotation(Dialog(group="Plot settings", colorSelector=true));

    output Real eigenvalues[size(A, 1), 2]
      "Eigenvalues of matrix A (Re: first column, Im: second column)";
  //  output Real evRe[size(A, 1)];
  //  output Real evIm[size(A, 1)];
  protected
    Integer markerStyle2=
      if markerStyle==MarkerStyles.Cross then MarkerStyle.Cross else
      if markerStyle==MarkerStyles.Circle then MarkerStyle.Circle else
      if markerStyle==MarkerStyles.Square then MarkerStyle.Square else
      if markerStyle==MarkerStyles.FilledSquare then MarkerStyle.FilledSquare else
      if markerStyle==MarkerStyles.TriangleDown then MarkerStyle.TriangleDown else
      if markerStyle==MarkerStyles.TriangleUp then MarkerStyle.TriangleUp else
      if markerStyle==MarkerStyles.Diamond then MarkerStyle.Diamond else MarkerStyle.Circle;
  //    if markerStyle==MarkerStyles.FilledCircle then MarkerStyle.FilledCircle else
    Boolean ok "True, if all calls are ok";

  algorithm
    // Calculate eigenvalues of input matrix A
    eigenvalues :=Modelica.Math.Matrices.eigenValues(A);
  //  evRe[:] :=eigenvalues[:, 1];
  //  evIm[:] :=eigenvalues[:, 2];

    // Plot real and imaginary part of eigenvalues
    if removePrevious then
      ok := removePlots();
      createPlot(id= 1,
        position=position,
        leftTitleType = 2,
        leftTitle = "Im",
        bottomTitleType = 2,
        bottomTitle = "Re",
        autoerase= false,
        grid=grid,
        heading=heading,
        legend=useLegend,
        erase= false,
        legendLocation = 2,
        legendHorizontal = false);
    end if;
    plotArray(
      x= eigenvalues[:, 1],
      y= eigenvalues[:, 2],
      legend=legend,
      color=markerColor,
      pattern = LinePattern.None,
      marker = markerStyle2);

  //   removePlots();
  // function createPlot "Create plot window"
  //   input Integer id := 0 "Window id";
  //   input Integer position[4] "Window Position";
  //   input String x := "time" "Independent variable";
  //   input String y[:] "Variables";
  //   input String heading := "" "Plot heading";
  //   input Real range[4] := {0.0, 1.0, 0.0, 1.0} "Range";
  //   input Boolean erase := true "Start with a fresh window";
  //   input Boolean autoscale := true "Autoscaling of y-axis";
  //   input Boolean autoerase := true "Erase previous when replotting";
  //   input Boolean autoreplot := true "Replot after simulation";
  //   input Boolean description := false "Include description in label";
  //   input Boolean grid := false "Add grid";
  //   input Boolean color := true "Deprecated. Replaced by colors, patterns, markers, and thicknesses.";
  //   input Boolean online := false "Online plotting";
  //   input Boolean legend := true "Variable legend";
  //   input Real timeWindow := 0.0 "Time window for online plotting";
  //   input String filename := "" "Result file to read data from";
  //   input Integer legendLocation := 1 "Where to place legend (1 above, 2 right, 3 below, 4-7 inside)";
  //   input Boolean legendHorizontal := true "Horizontal legend";
  //   input Boolean legendFrame := false "Draw frame around legend";
  //   input Boolean suppressMarker := false "Deprecated. Replaced by colors, patterns, markers, and thicknesses.";
  //   input Boolean logX := false "Logarithmic X scale";
  //   input Boolean logY := false "Logarithmic Y scale";
  //   input String legends[size(y, 1)] "Legends";
  //   input Integer subPlot := 1 "Sub plot number";
  //   input Boolean uniformScaling := false "Same vertical and horizontal axis increment";
  //   input Integer leftTitleType := 0 "Type of left axis title (0=none, 1=description, 2=custom)";
  //   input String leftTitle := "" "Custom left axis title";
  //   input Integer bottomTitleType := 0 "Type of bottom axis title (0=none, 1=description, 2=custom)";
  //   input String bottomTitle := "" "Custom bottom axis title";
  //   input Integer colors[size(y, 1), 3] "Line colors";
  //   input Integer patterns[size(y, 1)] "Line patterns, e.g., LinePattern.Solid";
  //   input Integer markers[size(y, 1)] "Line markers, e.g., MarkerStyle.Cross";
  //   input Real thicknesses[size(y, 1)] "Line thicknesses";
  //   output Integer _window;

    annotation (__Dymola_interactive=true);
  end plotEigenvalues;
end RootLocusOld;
