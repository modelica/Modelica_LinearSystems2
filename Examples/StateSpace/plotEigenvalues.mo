within Modelica_LinearSystems2.Examples.StateSpace;
function plotEigenvalues
  "Calculate eigenvalues of matrix A and plot root locii"

  input Real A[:,size(A, 1)] = [2,1,1;1,1,1;1,2,2] "Square matrix";
  input Boolean removePrevious=true
    "True, if all previous plots should be deleted";
  input String heading="Root locii" "Heading of plot";
  input String legend="Eigenvalues" "Legend of plot";
  input Integer markerStyle=1 "Style of marker"
    annotation (Dialog(group="Plot settings"),
      choices(
        choice=1 "Cross",
        choice=2 "Circle",
        choice=3 "Square",
        choice=4 "FilledSquare",
        choice=5 "TriangleDown",
        choice=6 "TriangleUp",
        choice=7 "Diamond"));
  input Integer markerColor[3]={0,0,255} "Color of marker";

  output Real eigenvalues[size(A, 1), 2]
    "Eigenvalues of matrix A (Re: first column, Im: second column)";
//  output Real evRe[size(A, 1)];
//  output Real evIm[size(A, 1)];
protected
  Integer markerStyle2=
    if markerStyle==1 then MarkerStyle.Cross else
    if markerStyle==2 then MarkerStyle.Circle else
    if markerStyle==3 then MarkerStyle.Square else
    if markerStyle==4 then MarkerStyle.FilledSquare else
    if markerStyle==5 then MarkerStyle.TriangleDown else
    if markerStyle==6 then MarkerStyle.TriangleUp else
    if markerStyle==7 then MarkerStyle.Diamond else MarkerStyle.Circle;
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
      leftTitleType=  2,
      leftTitle=  "Im",
      bottomTitleType=  2,
      bottomTitle=  "Re",
      autoerase= false,
      heading=heading,
      erase= false,
      legendLocation=  2,
      legendHorizontal=  false,
      position=  {35, 30, 665, 488});
  end if;
  plotArray(
    x= eigenvalues[:, 1],
    y= eigenvalues[:, 2],
    legend=legend,
    pattern=  LinePattern.None,
    marker=  markerStyle2);

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
//   input Boolean supressMarker := false "Deprecated. Replaced by colors, patterns, markers, and thicknesses.";
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

  annotation ();
end plotEigenvalues;
