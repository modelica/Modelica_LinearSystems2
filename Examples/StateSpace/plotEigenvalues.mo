within Modelica_LinearSystems2.Examples.StateSpace;
function plotEigenvalues
  "Linearize model, calculate eigenvalues and plot root locii"

  input Real A[:,size(A, 1)] = [2,1,1;1,1,1;1,2,2] "Square matrix";
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

  output Real eigenvalues[size(A, 1), 2]
    "Eigenvalues of matrix A (Re: first column, Im: second column)";
  output Real xxx[size(A, 1)];
  output Real yyy[size(A, 1)];
//   Real x[size(A, 1)];
//   Real y[size(A, 1)];
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
  eigenvalues :=Modelica.Math.Matrices.eigenValues(A);
  xxx[:] :=eigenvalues[:, 1];
  yyy[:] :=eigenvalues[:, 2];
  ok := removePlots();
  // createPlot(id=  1,
  //   filename=  "Unnamed.mat",
  //   autoerase=  false,
  //   position=  {35, 30, 665, 488});
  // plot(
  //   y=  {"x","y"},
  //   patterns=  {LinePattern.None,LinePattern.None},
  //   markers=  {MarkerStyle.Circle,MarkerStyle.Square});

  createPlot(id=  2,
    leftTitleType=  2,
    leftTitle=  "Im",
    bottomTitleType=  2,
    bottomTitle=  "Re",
    position=  {35, 30, 665, 488});
  plotArray(
    x= xxx,
    y=  yyy,
    legend="Root locii of ...",
    pattern=  LinePattern.None,
    marker=  markerStyle2);

//   removePlots();
//   createPlot(id=  2,
//     position=  {35, 30, 665, 488},
//     x=  "x",
//     y=  {"y"},
//     range=  {0.0, 1.0, 0.9, -0.1},
//     autoscale=  true,
//     autoerase=  true,
//     autoreplot=  true,
//     description=  false,
//     grid=  true,
//     color=  true,
//     online=  false,
//     filename=  "Unnamed.mat",
//     leftTitleType=  1,
//     bottomTitleType=  1,
//     colors=  {{0,0,255}},
//     patterns=  {LinePattern.None},
//     markers=  {MarkerStyle.Square});
    //MarkerStyle.None
    //MarkerStyle.Cross
    //MarkerStyle.Circle
    //MarkerStyle.Square
    //MarkerStyle.FilledSquare
    //MarkerStyle.TriangleDown
    //MarkerStyle.TriangleUp
    //MarkerStyle.Diamond

  annotation ();
end plotEigenvalues;
