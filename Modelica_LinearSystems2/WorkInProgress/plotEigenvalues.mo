within Modelica_LinearSystems2.WorkInProgress;
function plotEigenvalues
  "Calculate eigenvalues of matrix A and plot root locus"

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
  input MarkerStyle markerStyle=MarkerStyle.Cross "Style of marker"
    annotation (Dialog(group="Plot settings"));
  input Integer markerColor[3]={0,0,255} "Color of marker"
    annotation(Dialog(group="Plot settings", colorSelector=true));

  output Real eigenvalues[size(A, 1), 2]
    "Eigenvalues of matrix A (Re: first column, Im: second column)";
//  output Real evRe[size(A, 1)];
//  output Real evIm[size(A, 1)];
protected
  Integer id;
  Boolean ok "True, if all calls are ok";

algorithm
  // Calculate eigenvalues of input matrix A
  eigenvalues :=Modelica.Math.Matrices.eigenValues(A);
//  evRe[:] :=eigenvalues[:, 1];
//  evIm[:] :=eigenvalues[:, 2];

  // Plot real and imaginary part of eigenvalues
  if removePrevious then
    ok := DymolaCommands.Plot.removePlots();
    id := DymolaCommands.Plot.createPlot(id=1,
                    position=position,
                    y={""},
                    heading=heading,
                    erase=false,
                    autoscale=true,
                    autoerase=false,
                    grid=grid,
                    legend=useLegend,
                    legendLocation=2,
                    legendHorizontal=false,
                    legends={""},
                    leftTitleType = 2,
                    leftTitle = "Im",
                    bottomTitleType = 2,
                    bottomTitle = "Re",
                    colors=fill({-1, -1, -1}, 1),
                    patterns={LinePattern.Solid},
                    markers={MarkerStyle.None},
                    thicknesses={0.25},
                    axes={1},
                    displayUnits={""});
  end if;
  ok := DymolaCommands.Plot.plotArray(
    x= eigenvalues[:, 1],
    y= eigenvalues[:, 2],
    legend=legend,
    color=markerColor,
    pattern = LinePattern.None,
    marker = markerStyle,
    erase=false);

end plotEigenvalues;
