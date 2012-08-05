within Modelica_LinearSystems2.WorkInProgress;
package RootLocusOld
  function rootLocus
    "Plot root locus of nonlinear Modelica model by linearizing the model for variations of one model parameter"
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;
    import Modelica_LinearSystems2.Utilities.Types.MarkerStyles;

    input String modelName "Name of the Modelica model"
      annotation(Dialog(__Dymola_translatedModel(caption="Model to be linearized for the root locus")));

    input Modelica.SIunits.Time t_linearize= 0
      "Simulate until t_linearize and then linearize" annotation(Dialog);

  /*
  input Modelica.Mechanics.MultiBody.Interfaces.partialColorMap colorMap=
      Modelica.Mechanics.MultiBody.Visualizers.Colors.ColorMaps.jet
    "Color map function for parameter variations" annotation(__Dymola_choicesAllMatching=true);
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
    Boolean grid=true "True, if grid is shown" annotation(Dialog,  choices(__Dymola_checkBox=true));

    /* group "Axes" (Axes properties) */
    String xLabel="Real part of eigenvalues"
      "String displayed at horizontal axis"                                        annotation(Dialog(group="Axes"));
    String yLabel="Imaginary part of eigenvalues"
      "String displayed at vertical axis"                                             annotation(Dialog(group="Axes"));
    Boolean logX = false "True, if logarithmic scale of x-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
    Boolean logY = false "True, if logarithmic scale of y-axis" annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));
    Boolean uniformScaling = false
      "True, if same vertical and horizontal axis increment"
        annotation(Dialog(group="Axes"),choices(__Dymola_checkBox=true));

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
    Real nxMat[1,1]=readMatrix(fileName2, "nx", 1, 1);
    Integer ABCDsizes[2]=readMatrixSize(fileName2, "ABCD");
    Integer nx=integer(nxMat[1, 1]);
    Integer nu=ABCDsizes[2] - nx;
    Integer ny=ABCDsizes[1] - nx;
    Real ABCD[nx + ny,nx + nu]=readMatrix(fileName2, "ABCD", nx + ny, nx + nu);
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
end RootLocusOld;
