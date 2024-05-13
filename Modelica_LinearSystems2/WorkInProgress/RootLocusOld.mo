within Modelica_LinearSystems2.WorkInProgress;
package RootLocusOld
  function rootLocus
    "Plot root locus of nonlinear Modelica model by linearizing the model for variations of one model parameter"
    import Modelica_LinearSystems2;
    import Modelica_LinearSystems2.StateSpace;

    input String modelName "Name of the Modelica model"
      annotation(Dialog(__Dymola_translatedModel(caption="Model to be linearized for the root locus")));

    input Modelica.Units.SI.Time t_linearize=0
      "Simulate until t_linearize and then linearize" annotation (Dialog);

  /*
  input Modelica.Mechanics.MultiBody.Interfaces.partialColorMap colorMap=
      Modelica.Mechanics.MultiBody.Visualizers.Colors.ColorMaps.jet
    "Color map function for parameter variations" annotation(choicesAllMatching=true);
*/

    input
      Modelica_LinearSystems2.Records.ParameterVariation modelParam[:]
      "Model parameter to be varied";

    input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
        Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
      "Simulation options it t_linearize > 0";

    input
      Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram diagram
      "Diagram properties of the root locus";
    input Modelica_LinearSystems2.Utilities.Plot.Records.Device device=
      Modelica_LinearSystems2.Utilities.Plot.Records.Device()
      "Properties of device where figure is shown" annotation(Dialog);

  protected
    Boolean ok "True, if call is ok";
    Real dp "Step of parameter equidistand grid";
    Real parValue "Value of parameter in a loop";
    Integer nVar = max(modelParam[1].nPoints,2)
      "Number of parameter variations (at least 2 are necessary)";
    Real color[nVar,3] = Modelica.Mechanics.MultiBody.Visualizers.Colors.ColorMaps.jet(nVar)
      "Color of markers in a loop";
    Real mmToPixel= device.windowResolution/25.4;
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
    position := {device.xTopLeft,
                 device.yTopLeft,
                 device.diagramWidth,
                 diagram.heightRatio*device.diagramWidth}*mmToPixel;
    id:= DymolaCommands.Plot.createPlot(id=-1,
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
    ok:=DymolaCommands.SimulatorAPI.translateModel(modelName);
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
      ok := DymolaCommands.Plot.plotArray(eigenValues[:,1],
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
Plot.<strong>rootLocus</strong>(modelName, t_linearize, modelParam, simulationSetup, diagram)
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
Utilities.Plot.<strong>rootLocus</strong>(
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
<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/>
</div>
</html>"));
  end rootLocus;

  function rootLocusOfDriveOld
    "Plot the root locus of a drive with varying load"
  algorithm
    Modelica_LinearSystems2.WorkInProgress.RootLocusOld.rootLocus(
      "Modelica.Mechanics.Rotational.Examples.First",
      modelParam={
        Modelica_LinearSystems2.Records.ParameterVariation(Name="Jload", Min=1, Max=20, nPoints=30)},
      diagram=Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram(
        heading="Root loci of simple drive train"));
      annotation(__Dymola_interactive=true, Documentation(info="<html>
<p>
This function plots the root locus of model
<a href=\"modelica://Modelica.Mechanics.Rotational.Examples.First\">Rotational.Examples.First</a>
over the load inertia <strong>Jload</strong>:
</p>

<div>
<img src=\"modelica://Modelica_LinearSystems2/Resources/Images/Utilities/rootLocusOfDrive.png\">
</div>
</html>"));
  end rootLocusOfDriveOld;

  function linearize2
    "Linearize a model after simulation up to a given time and return only the A matrix"
    import Simulator = DymolaCommands.SimulatorAPI;

    input String modelName "Name of the Modelica model" annotation(Dialog(__Dymola_translatedModel));
    input Modelica.Units.SI.Time t_linearize=0
      "Simulate until t_linearize and then linearize" annotation (Dialog);
    input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
        Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
      "Simulation options it t_linearize > 0" annotation(Dialog);

  protected
    String fileName="dslin";
    String fileName2=fileName+".mat";

    // Simulate until t_linearize and then linearize at this time instant
    Boolean OK1 = if t_linearize <= 0.0 then true else
                     Simulator.simulateModel(problem=modelName, startTime=0, stopTime=t_linearize,
                                   method=simulationSetup.method,
                                   tolerance=simulationSetup.tolerance,
                                   fixedstepsize=simulationSetup.fixedStepSize);
    Boolean OK2 = if t_linearize <= 0.0 then true else Simulator.importInitial("dsfinal.txt");
    Boolean OK3 = Simulator.linearizeModel(problem=modelName, resultFile=fileName, startTime=t_linearize, stopTime=t_linearize);

    // Read linear system from file
    Integer xuy[3]=Modelica_LinearSystems2.Utilities.Streams.readSystemDimension(
      fileName2, "ABCD");
    Integer nx = xuy[1];
    Integer nu = xuy[2];
    Integer ny = xuy[3];
    Real ABCD[nx + ny,nx + nu]=Modelica.Utilities.Streams.readRealMatrix(fileName2, "ABCD", nx + ny, nx + nu);
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
