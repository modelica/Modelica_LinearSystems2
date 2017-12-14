within Modelica_LinearSystems2.WorkInProgress.Tests.DiscreteSystems;
function Operations

  import Modelica;
  import Modelica_LinearSystems2;
  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Modelica_LinearSystems2.TransferFunction;
  import Modelica_LinearSystems2.DiscreteStateSpace;
  import Modelica_LinearSystems2.DiscreteZerosAndPoles;
  import Modelica_LinearSystems2.DiscreteTransferFunction;
  import Modelica.Math.Vectors;

  input Boolean doPlot=true;
  input Modelica.SIunits.Time Ts=0.1;

protected
  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact;

  ZerosAndPoles zp1=ZerosAndPoles(
      n1={2},
      n2=[2,2],
      d1={1,1.5},
      d2=[1,1; 2,3]);
  ZerosAndPoles zp2=zp1 + zp1;
  DiscreteZerosAndPoles dzp1=DiscreteZerosAndPoles(zp1, Ts, method);
  DiscreteStateSpace dss1=DiscreteStateSpace(dzp1);
  DiscreteTransferFunction dtf1=
      DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp1);

  DiscreteStateSpace dss2=dss1 + dss1;
  DiscreteTransferFunction dtf2=dtf1 + dtf1;
  DiscreteZerosAndPoles dzp2=dzp1 + dzp1;

  DiscreteZerosAndPoles dzp3=2*dzp1 - dzp1;
  DiscreteTransferFunction dtf3=2*dtf1 - dtf1;
  DiscreteStateSpace dss3=2*dss1-dss1;

  DiscreteZerosAndPoles dzp4=dzp3/(1/dzp1);
  DiscreteTransferFunction dtf4=dtf3/(1/dtf1);

  Real tSpan=10;
  Real y[:,1,1];
  Real y0[:];
  Real y1[:];
  Real y2[:];
  Real t[:];
  Real delta;
  Boolean ok=true;

algorithm
  if doPlot then
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss2, tSpan=tSpan);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp2, tSpan=tSpan);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.step(dtf2, tSpan=tSpan);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss3, tSpan=tSpan);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp3, tSpan=tSpan);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.step(dtf3, tSpan=tSpan);
  end if;

// ###############  '+'  #################
  y := DiscreteZerosAndPoles.Analysis.stepResponse(dzp=dzp2, tSpan=tSpan);
  y0 := y[:, 1, 1];

  y := DiscreteStateSpace.Analysis.stepResponse(dss=dss2, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/Vectors.norm(y0)/size(y0, 1);
  Modelica.Utilities.Streams.print("delta_1 = " + String(delta));
  ok := delta < 1e-6;
  assert(ok, "dzp2 or dss2 failed");

  y := DiscreteTransferFunction.Analysis.stepResponse(dtf=dtf2, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/Vectors.norm(y0)/size(y0, 1);
  Modelica.Utilities.Streams.print("delta_2 = " + String(delta));
  ok := delta < 1e-4;
  assert(ok, "dzp2 or dtf2 failed");

// ###############  '-'  #################
  y := DiscreteZerosAndPoles.Analysis.stepResponse(dzp=dzp3, tSpan=tSpan);
  y0 := y[:, 1, 1];

  y := DiscreteStateSpace.Analysis.stepResponse(dss=dss3, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/Vectors.norm(y0)/size(y0, 1);
  Modelica.Utilities.Streams.print("delta_3 = " + String(delta));
  ok := delta < 1e-6;
  assert(ok, "dzp3 or dss3 failed");

  (y,t) := DiscreteTransferFunction.Analysis.stepResponse(dtf=dtf3, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y0, 1);
  Modelica.Utilities.Streams.print("delta_4 = " + String(delta));
  ok := delta < 1e-4;
//  assert(ok, "dzp3 or dtf3 failed");

  Modelica_LinearSystems2.Utilities.Plot.diagram(
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
        curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=t,
          y=y0,
          legend="y0"),
      Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
          x=t,
          y=y1,
          legend="y1")}));
// ###############  '/'  #################
  (y,t) := DiscreteZerosAndPoles.Analysis.stepResponse(dzp=dzp4, tSpan=6);
  y0 := y[:, 1, 1];

  y := DiscreteTransferFunction.Analysis.stepResponse(dtf=dtf4, tSpan=6);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y0, 1);
  Modelica.Utilities.Streams.print("delta_5 = " + String(delta));
  ok := delta < 1e-4;
//  assert(ok, "dzp4 or dtf4 failed");

  Modelica_LinearSystems2.Utilities.Plot.diagram(
    Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(curve={
    Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
    x=t,
    y=y0,
    legend="y0"),Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
    x=t,
    y=y1,
    legend="y1")}));

  Modelica.Utilities.Streams.print("\n OK !!!! ");

  annotation(__Dymola_interactive=true);
end Operations;
