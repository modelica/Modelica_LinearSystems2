within Modelica_LinearSystems2.Utilities.Plot;
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

  input Modelica_LinearSystems2.Records.ParameterVariation modelParam[:]
    "Model parameter to be varied";

  input Modelica_LinearSystems2.Records.SimulationOptionsForLinearization simulationSetup=
      Modelica_LinearSystems2.Records.SimulationOptionsForLinearization()
    "Simulation options it t_linearize > 0";

  input Modelica_LinearSystems2.Utilities.Plot.Records.RootLocusDiagram diagram
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

  // Loop over the selected parameter range
  for i in 1:nVar loop
     // Linearize
     parValue := modelParam[1].Min+(i-1)*dp;
     ok :=SetVariable(modelParam[1].Name, parValue);
     A :=Modelica_LinearSystems2.Utilities.Import.linearize2(
             modelName,
             t_linearize=t_linearize,
             simulationSetup=simulationSetup);

     // compute eigen values
     eigenValues :=Modelica.Math.Matrices.eigenValues(A);

     // Plot root loci
     ok :=plotArray(eigenValues[:,1],
                    eigenValues[:,2],
                    legend=String(parValue,significantDigits=3),
                    color=integer(color[i, :]),
                    pattern=LinePattern.None,
                    marker=MarkerStyle.Square,
                    id=id);
  end for;
  annotation (__Dymola_interactive=true, Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>plotRootLoci(modelName, modelParam, simulationOptions, position, useLegend, grid, markerStyle, markerColorMin, markerColorMax)</pre></blockquote>

<h4>Description</h4>
<p>Calling this function with default setup the following root loci plot will be generated.</p>
<p><img src=\"modelica://Modelica_LinearSystems2/Resources/Images/rootLociiDefaultSetup.png\"/></p>
</html>"));
end rootLocus;
