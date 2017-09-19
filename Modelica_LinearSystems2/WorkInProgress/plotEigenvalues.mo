within Modelica_LinearSystems2.WorkInProgress;
function plotEigenvalues
  "Calculate eigenvalues of matrix A and plot root locus"
  import Modelica_LinearSystems2.WorkInProgress.RootLocusOld.Types.MarkerStyles;

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
  Integer id;
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
    id := createPlot(id= 1,
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
  ok := plotArray(
    x= eigenvalues[:, 1],
    y= eigenvalues[:, 2],
    legend=legend,
    color=markerColor,
    pattern = LinePattern.None,
    marker = markerStyle2,
    erase=false);

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
