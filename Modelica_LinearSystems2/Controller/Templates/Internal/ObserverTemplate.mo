within Modelica_LinearSystems2.Controller.Templates.Internal;
block ObserverTemplate
  "Template of a Luenberger observer for state space systems"
  extends Interfaces.PartialSampledBlock;

  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.Internal.StateSpace2;

  parameter Boolean matrixOnFile=false
    "True, if matrix should be read from file";
  parameter String fileName=Modelica_LinearSystems2.DataDir + "observer.mat"
    "Name of the state space system data file"
    annotation(Dialog(loadSelector(filter="MAT files (*.mat);; All files (*.*)",
      caption="observer data file"),enable = matrixOnFile));
  parameter String systemName="stateSpace" "Name of state space system" annotation(Dialog(enable = matrixOnFile));
  parameter String observerMatrixName="L" "Name of matrix" annotation(Dialog(enable = matrixOnFile));

  parameter StateSpace2 plantModelSystem=StateSpace2(A=[0],B=[1],C=[1],D=[0])
    "Plant state space system" annotation(Dialog(enable = not matrixOnFile));
  parameter Real L[:,:]=[1] "Observer feedback matrix" annotation(Dialog(enable = not matrixOnFile));

protected
  parameter Integer mn[2]=if matrixOnFile then
    Modelica.Utilities.Streams.readMatrixSize(fileName, observerMatrixName) else
    size(L);
  parameter Integer m=mn[1];
  parameter Integer n=mn[2];

  parameter Real L2[:,:]=if matrixOnFile then
    Modelica.Utilities.Streams.readRealMatrix(fileName, observerMatrixName, m, n) else L;
  parameter StateSpace2 plantModelSystem2=
    if matrixOnFile then StateSpace2.Import.fromFile(
    fileName, systemName) else plantModelSystem;
  parameter Real C[:,:]=plantModelSystem2.C;

public
  parameter Boolean withDelay = Modelica_LinearSystems2.Math.Matrices.Internal.haveZeroRow(observerStateSpace.system.A)
    "True, if a unit delay should be considered";

  final parameter Integer nout=size(L2,2)
    "Dimension of vector of measured output (y)";
  parameter Real x_start[size(plantModelSystem2.A,1)]=zeros(size(plantModelSystem2.A,1))
    "Initial or guess values of states" annotation(Dialog(tab="Advanced options"));
//  parameter Real y_start[size(plantModelSystem2.C,1)]=zeros(size(plantModelSystem2.C,1)) "Initial values of outputs (remaining states are in steady state if possible)" annotation 6;
 //    initType= if init==Types.Init.InitialState then Types.Init.InitialState else  Types.Init.NoInit
// y_start-plantModelSystem2.C*L2*observerStateSpace.y_start

  Modelica.Blocks.Interfaces.RealInput u[size(plantModelSystem2.B, 2)]
    annotation (Placement(transformation(extent={{-140,40},{-100,80}})));
  Modelica.Blocks.Interfaces.RealInput y[nout] annotation (Placement(transformation(extent={{-140,-80},{-100,-40}})));
  Modelica.Blocks.Interfaces.RealOutput x_estimated[size(observerStateSpace.system.A,1)]
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));
  Modelica.Blocks.Routing.Multiplex2 multiplex2_1(n1=size(plantModelSystem2.B,2), n2=nout)
    annotation (Placement(transformation(extent={{-60,-10},{-40,10}})));
  Controller.StateSpace observerStateSpace(
    x_start= x_start,
    initType=Types.InitWithGlobalDefault.NoInit,
    system(
      A=plantModelSystem2.A - L2*plantModelSystem2.C,
      B=[plantModelSystem2.B,L2],
      C=identity(size(plantModelSystem2.A, 1)),
      D=zeros(size(plantModelSystem2.A, 1), size(plantModelSystem2.B, 2)+nout)),
    withDelay=withDelay)
    annotation (Placement(transformation(extent={{-20,-20},{20,20}})));
initial equation
  //  if continuous then
    if init == Types.Init.InitialState then
      x_estimated = x_start-L2*plantModelSystem2.C*x_start;
    elseif init == Types.Init.SteadyState then
      der(x_estimated) = zeros(size(plantModelSystem2.A,1));
    elseif init == Types.Init.InitialOutput then
      x_estimated=Modelica.Math.Matrices.inv(transpose(C)*C)*transpose(C)*x_start-L2*plantModelSystem2.C*x_start;
      der(x_estimated[size(plantModelSystem2.B,2)+1:size(plantModelSystem2.A,1)]) = zeros(size(plantModelSystem2.A,1)-size(plantModelSystem2.B,2));
    end if;
//  end if;

equation
  connect(multiplex2_1.u1, u) annotation (Line(
      points={{-62,6},{-72,6},{-72,60},{-120,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiplex2_1.u2, y) annotation (Line(
      points={{-62,-6},{-72,-6},{-72,-60},{-120,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(x_estimated, observerStateSpace.y) annotation (Line(
      points={{110,0},{22,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(observerStateSpace.u, multiplex2_1.y) annotation (Line(
      points={{-24,0},{-39,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics={
        Text(
          extent={{-138,46},{-78,30}},
          textColor={0,0,255},
          textString="(control input)"),
        Text(
          extent={{-140,-72},{-80,-88}},
          textColor={0,0,255},
          textString="(measured output)"),
        Text(
          extent={{90,28},{150,12}},
          textColor={0,0,255},
          textString="(estimated state)"),
        Text(
          extent={{-60,-20},{60,-70}},
          textColor={0,0,255},
          textString="der(x) = (A-LC)x + Bu + Ly"),
        Text(
          extent={{-50,-52},{-28,-64}},
          textColor={0,0,255},
          textString="y = x")}),
    Icon(graphics={
        Polygon(
          points={{-78,92},{-86,70},{-70,70},{-78,92}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-78,80},{-78,-88}}, color={192,192,192}),
        Line(points={{-88,-78},{84,-78}}, color={192,192,192}),
        Polygon(
          points={{92,-78},{70,-70},{70,-86},{92,-78}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-78,28},{-56,38},{-30,30},{-4,22},{30,36},{46,34},{66,22},
              {76,10},{84,6}},
          color={0,0,255},
          smooth=Smooth.None),
        Line(
          points={{-78,-30},{-54,10},{-30,18},{-4,16},{32,30},{48,32},{66,20},
              {76,10},{84,6}},
          color={255,0,0},
          smooth=Smooth.None),
        Text(
          extent={{-48,4},{-24,-18}},
          textColor={255,0,0},
          textString="x"),
        Text(
          extent={{-46,10},{-24,-10}},
          textColor={255,0,0},
          textString="^"),
        Text(
          extent={{-76,-30},{76,-76}},
          textColor={0,0,255},
          textString="state estimation"),
        Text(
          extent={{-58,58},{-34,36}},
          textColor={0,0,255},
          textString="x"),
        Text(
          extent={{-128,94},{-108,76}},
          textColor={0,0,255},
          textString="u"),
        Text(
          extent={{-128,-28},{-108,-46}},
          textColor={0,0,255},
          textString="y")}));
end ObserverTemplate;
