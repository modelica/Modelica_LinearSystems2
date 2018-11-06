within Modelica_LinearSystems2.WorkInProgress.Tests.DiscreteSystems;
function Conversions

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
  input String modelName="Modelica_LinearSystems2.Utilities.Plants.DoublePendulum";
  input Modelica.SIunits.Time Ts=0.1;

protected
  Modelica_LinearSystems2.Utilities.Types.Method method=Modelica_LinearSystems2.Utilities.Types.Method.StepExact;

  ZerosAndPoles zp1=ZerosAndPoles(
      n1={2},
      n2=[2,2],
      d1={1,1.5},
      d2=[1,1; 2,3]);
  TransferFunction tf1=ZerosAndPoles.Conversion.toTransferFunction(zp1);
  StateSpace ss1=StateSpace(zp1);
  DiscreteZerosAndPoles dzp1=DiscreteZerosAndPoles(zp1, Ts, method);
  DiscreteStateSpace dss1=DiscreteStateSpace(dzp1);
  DiscreteTransferFunction dtf1=DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp1);
  DiscreteStateSpace dss2=DiscreteStateSpace(ss1, Ts, method);
  DiscreteTransferFunction dtf2=DiscreteTransferFunction(tf1, Ts, method);
  DiscreteZerosAndPoles dzp2=DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss2);
  DiscreteZerosAndPoles dzp3=DiscreteZerosAndPoles(dtf2);
  DiscreteTransferFunction dtf3=DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss2);
  DiscreteStateSpace dss3=DiscreteStateSpace(dtf2);

  DiscreteStateSpace dss4=Modelica_LinearSystems2.DiscreteStateSpace.Import.fromModel(modelName=modelName,Ts=Ts);
  DiscreteZerosAndPoles dzp4[6,1]=Modelica_LinearSystems2.DiscreteZerosAndPoles.Import.fromModel(modelName=modelName,Ts=Ts);
  DiscreteTransferFunction dtf4[6,1]=Modelica_LinearSystems2.DiscreteTransferFunction.Import.fromModel(modelName=modelName,Ts=Ts);

  Real tSpan=10;
  Real y[:,1,1];
  Real y0[:];
  Real y1[:];
  Real y2[:];
  Real delta;
  Boolean ok=true;

algorithm
  if doPlot then
    Modelica_LinearSystems2.ZerosAndPoles.Plot.step(zp1);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss1);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss2);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss3);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp1);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp2);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp3);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.step(dtf1);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.step(dtf2);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.step(dtf3);
    Modelica_LinearSystems2.ZerosAndPoles.Plot.ramp(zp1);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp(dss1);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp(dss2);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.ramp(dss3);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp(dzp1);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp(dzp2);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.ramp(dzp3);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.ramp(dtf1);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.ramp(dtf2);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.ramp(dtf3);

//    Modelica_LinearSystems2.DiscreteStateSpace.Plot.step(dss4);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.step(dzp4[1,1]);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.step(dtf4[1,1]);
  end if;

// ###############   step - sesponses  #################
  y := DiscreteZerosAndPoles.Analysis.stepResponse(dzp=dzp1, tSpan=tSpan);
  y0 := y[:,1,1];
  y := DiscreteZerosAndPoles.Analysis.stepResponse(dzp=dzp2, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_1 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dzp1 or dzp2 failed");
  y := DiscreteZerosAndPoles.Analysis.stepResponse(dzp=dzp3, tSpan=tSpan);
  y2 := y[:, 1, 1];
  delta := Vectors.norm(y1 - y2)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_2 = "+String(delta));
  ok := delta<1e-9;
  assert(ok,"dzp1 or dzp3 failed");

  y := DiscreteStateSpace.Analysis.stepResponse(dss=dss1, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_3 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dzp1 or dss1 failed");
  y := DiscreteStateSpace.Analysis.stepResponse(dss=dss2, tSpan=tSpan);
  y2 := y[:, 1, 1];
  delta := Vectors.norm(y1 - y2)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_4 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dss1 or dss2 failed");
  y := DiscreteStateSpace.Analysis.stepResponse(dss=dss3, tSpan=tSpan);
  y2 := y[:, 1, 1];
  delta := Vectors.norm(y1 - y2)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_5 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dss1 or dss3 failed");

  y := DiscreteTransferFunction.Analysis.stepResponse(dtf=dtf1, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_6 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dzp1 or dtf1 failed");
  y := DiscreteTransferFunction.Analysis.stepResponse(dtf=dtf2, tSpan=tSpan);
  y2 := y[:, 1, 1];
  delta := Vectors.norm(y1 - y2)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_7 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dtf1 or dtf2 failed");
  y := DiscreteTransferFunction.Analysis.stepResponse(dtf=dtf3, tSpan=tSpan);
  y2 := y[:, 1, 1];
  delta := Vectors.norm(y1 - y2)/size(y1,1);
  Modelica.Utilities.Streams.print("delta_8 = "+String(delta));
  ok := delta<1e-10;
  assert(ok,"dtf2 or dtf3 failed");

 // ########  ramp responses  ########
  y := DiscreteZerosAndPoles.Analysis.rampResponse(dzp=dzp1, tSpan=tSpan);
  y0 := y[:,1,1];

  y := DiscreteStateSpace.Analysis.rampResponse(dss=dss1, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1, 1);
  Modelica.Utilities.Streams.print("delta_9 = " + String(delta));
  ok := delta < 1e-10;
  assert(ok, "dzp1 or dss1 failed");

  y := DiscreteTransferFunction.Analysis.rampResponse(dtf=dtf1, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1, 1);
  Modelica.Utilities.Streams.print("delta_10 = " + String(delta));
  ok := delta < 1e-10;
  assert(ok, "dzp1 or dtf1 failed");

//##############   impulse check   ###################
  method :=Modelica_LinearSystems2.Utilities.Types.Method.ImpulseExact;

  dzp1 := DiscreteZerosAndPoles(zp1, Ts, method);
  dss1 := DiscreteStateSpace(dzp1);
  dtf1 := DiscreteZerosAndPoles.Conversion.toDiscreteTransferFunction(dzp1);
  dss2 := DiscreteStateSpace(ss1, Ts, method);
  dtf2 := DiscreteTransferFunction(tf1, Ts, method);
  dzp2 := DiscreteStateSpace.Conversion.toDiscreteZerosAndPoles(dss2);
  dzp3 := DiscreteZerosAndPoles(dtf2);
  dtf3 := DiscreteStateSpace.Conversion.toDiscreteTransferFunction(dss2);
  dss3 := DiscreteStateSpace(dtf2);

  y := DiscreteZerosAndPoles.Analysis.impulseResponse(dzp=dzp1, tSpan=tSpan);
  y0 := y[:, 1, 1];

  y := DiscreteStateSpace.Analysis.impulseResponse(dss=dss1, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1, 1);
  Modelica.Utilities.Streams.print("delta_11 = " + String(delta));
  ok := delta < 1e-10;
  assert(ok, "dzp1 or dss1 failed");

  y := DiscreteTransferFunction.Analysis.impulseResponse(dtf=dtf1, tSpan=tSpan);
  y1 := y[:, 1, 1];
  delta := Vectors.norm(y0 - y1)/size(y1, 1);
  Modelica.Utilities.Streams.print("delta_12 = " + String(delta));
  ok := delta < 1e-10;
  assert(ok, "dzp1 or dtf1 failed");

  if doPlot then
    Modelica_LinearSystems2.ZerosAndPoles.Plot.impulse(zp1);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse(dss1);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse(dss2);
    Modelica_LinearSystems2.DiscreteStateSpace.Plot.impulse(dss3);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse(dzp1);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse(dzp2);
    Modelica_LinearSystems2.DiscreteZerosAndPoles.Plot.impulse(dzp3);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.impulse(dtf1);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.impulse(dtf2);
    Modelica_LinearSystems2.DiscreteTransferFunction.Plot.impulse(dtf3);
  end if;

  DiscreteStateSpace.Plot.bodeSISO(dss2);
  DiscreteZerosAndPoles.Plot.bode(dzp2);
  DiscreteTransferFunction.Plot.bode(dtf2);

  Modelica.Utilities.Streams.print("ok = "+String(ok));

end Conversions;
