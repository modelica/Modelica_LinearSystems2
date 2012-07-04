within Modelica_LinearSystems2.WorkInProgress.Controller.Examples;
model ZerosAndPolesBlock
  extends Modelica.Icons.Example;

  import Modelica_LinearSystems2.ZerosAndPoles;

  parameter ZerosAndPoles zp=ZerosAndPoles(n1=fill(0,0), d1={0,0}, n2=fill(0,0,2), d2=fill(0,0,2));

//  parameter TransferFunction tf=ZerosAndPoles.Conversion.toTransferFunction(zp);

  Modelica_LinearSystems2.Controller.ZerosAndPoles zerosAndPoles(
    system(
    k=4096,
      n1 = {0, -3.0, -1, -1, -1, 1, 1, 1, 1.5, 1.5, 2, 2, 2, 2, 3.0, 4.0},
      n2 = [-1, 2.0;
1, 1.0;
1, 1.0;
1, 2.0],
      d1 =  {1, 1, 1, 1, 1, 1, 1, 1, 1.1, 1.1, 2, 2, 2, 2, 2, 2, 2.5, 2.6, 3.0, 3.3, 3.5},
      d2 =  [-1.0, 1.5;
0, 1;
0, 2.0;
1.0, 3.0;
1.0, 3.0;
2.0, 2.0;
2.0, 2.0;
2.0, 2.0;
2, 3.0]),
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption,
    initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.NoInit)
    annotation (Placement(transformation(extent={{-20,-40},{0,-20}})));

  Modelica_LinearSystems2.Controller.TransferFunction transferFunction(
    system(n = {4096.0, 69632.0, 484352.0, 1660928.0, 1910784.0, -7356416.0, -40869888.0,
  -102811648.0, -165557248.0, -172556288.0, -69864447.9999999, 159680512.0,
  464871424.0, 681885696.0, 565237760.0, 24422399.9999997, -634348544.0,
  -892178431.999999, -558878719.999999, 4603904.00000018, 318455808.0,
  281346048.0, 118554624.0, 21233664.0, 0}, d = {1, 46.1, 1032.82, 15010.0, 159357.8593, 1319621.81465, 8886051.5153,
  50084410.513025, 241347604.02505, 1010571342.0867, 3723992761.10015,
  12201626614.2633, 35845189862.7856, 95073444683.0585, 228990221232.589,
  503283393417.459, 1013474627929.38, 1876226721140.07, 3202024690492.75,
  5048540798554.17, 7365220535364.67, 9951666633258.59, 12457405523701.7,
  14441507868912.9, 15486650166621.3, 15332613929249.6, 13975117865117.0,
  11681867052184.8, 8911174949396.56, 6164516162114.78, 3836912153456.54,
  2127539809735.55, 1037814234428.29, 438182671104.399, 156735759450.941,
  46122743997.6048, 10703484314.0928, 1833242319.5328, 205601583.0528,
  11302042.752}),
    blockType=Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.UseSampleClockOption,
    initType=Modelica_LinearSystems2.Controller.Types.InitWithGlobalDefault.NoInit)
    annotation (Placement(transformation(extent={{-20,0},{0,20}})));

  Modelica.Blocks.Sources.Step step(startTime=1)
    annotation (Placement(transformation(extent={{-60,-20},{-40,0}})));
  inner Modelica_LinearSystems2.Controller.SampleClock sampleClock(
    methodType=Modelica_LinearSystems2.Types.Method.Trapezoidal,
    sampleTime=0.01,
    initType=Modelica_LinearSystems2.Controller.Types.Init.InitialOutput,
    blockType=Modelica_LinearSystems2.Controller.Types.BlockType.Discrete)
    annotation (Placement(transformation(extent={{60,60},{80,80}})));

equation
  connect(zerosAndPoles.u, step.y) annotation (Line(
      points={{-22,-30},{-32,-30},{-32,-10},{-39,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(transferFunction.u, step.y) annotation (Line(
      points={{-22,10},{-32,10},{-32,-10},{-39,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    experiment(StopTime=5),
    experimentSetupOutput);
end ZerosAndPolesBlock;
