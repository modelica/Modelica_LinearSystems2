within Modelica_LinearSystems2.Controller;
block Noise
  "Block that generates a uniform distributed noise output signal at sample instants if discrete and y=0 if continuous"
  extends Interfaces.PartialBlockIcon;
  parameter Real y_min if not continuous "Lower limit of noise band";
  parameter Real y_max if not continuous "Upper limit of noise band";
  parameter Integer firstSeed[3](
    each min=0,
    each max=255) = {23,87,187} if                                          not continuous
    "Integer[3] defining random sequence; required element range: 0..255";
  parameter Types.BlockTypeWithGlobalDefault blockType=Types.BlockTypeWithGlobalDefault.UseSampleClockOption
    "Type of block (Continuous/Discrete)" 
    annotation(Evaluate=true, Hide=true);
  final parameter Boolean continuous=blockType == Types.BlockTypeWithGlobalDefault.Continuous
       or blockType == Types.BlockTypeWithGlobalDefault.UseSampleClockOption
       and sampleClock.blockType == Types.BlockType.Continuous
    "= true, if continuous block, otherwise discrete block";
  parameter Integer sampleFactor(min=1) = 1 if 
                                             not continuous
    "Ts=sampleClock.sampleTime*sampleFactor" 
     annotation (Dialog(enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockTypeWithGlobalDefault.Continuous));
  Modelica.Blocks.Interfaces.RealOutput y "Discrete output signal of block" 
    annotation (extent=[100, -10; 120, 10]);

  annotation (
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Polygon(
          points={{-81,90},{-89,68},{-73,68},{-81,90}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-81,78},{-81,-90}}, color={192,192,192}),
        Line(points={{-90,-23},{82,-23}}, color={192,192,192}),
        Polygon(
          points={{91,-22},{69,-14},{69,-30},{91,-22}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-24,95},{94,56}},
          lineColor={192,192,192},
          textString="noise"),
        Line(points={{-35,13},{-35,-47},{-25,-47},{-25,-29},{-15,-29},{-15,-57},
              {-5,-57},{-5,25},{1,25},{1,39},{7,39},{7,-17},{17,-17},{17,-5},{
              23,-5},{23,-35},{33,-35},{33,37},{43,37},{43,3},{51,3},{51,-63},{
              61,-63}}, color={0,0,127}),
        Line(points={{-81,-29},{-67,-29},{-67,-13},{-59,-13},{-59,-61},{-51,-61},
              {-51,-39},{-43,-39},{-43,45},{-35,45},{-35,13}}, color={0,0,127}), 

        Line(points={{-90,-70},{84,-70}}, color={255,0,0}),
        Line(points={{-89,50},{85,50}}, color={255,0,0}),
        Text(
          extent={{-145,-105},{144,-133}},
          lineColor={0,0,0},
          textString="[%y_min .. %y_max]")}),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics),
    Documentation(info="<html>
<p>
If <b>discrete</b> block, the output y is sampled according to sample time
sampleClock.sampleTime * sampleFactor, where sampleClock.sampleTime
is defined globally in the outer component sampleClock and
sampleFactor is an Integer parameter of component Noise.
At every sample time, a random output signal y in the range y_min .. y_max
is generated, where y_min and y_max are parameters. A typical
noise signal is shown in the next figure:
</p>
<p align=\"center\">
<img src=\"../Extras/Images/Noise1.png\">
</p>
<p>
The Integer[3] parameter vector <b>firstSeed</b> is used to initialize the
basic random number generator. The 3 elements of firstSeed need
to be in the range [0, 255]. The use of the same seed vector 
will lead to the same sequence of numbers when these are computed serially.
This is usually not desired. Therefore, for every usage of block
<b>Noise</b> a different firstSeed should be defined.
</p>
<p>
If <b>continuous</b> block, the output y = 0.0, i.e., no noise signal
is generated. The reason is that the noise can only reasonably be
used in a simulation if it is a discrete signal, i.e., changes
its value only at sample instants. Since a continous block is usually
used to speed up the simulation, the noise should also be turned
off because it will otherwise significantly limit the
maximum step size of the integrator.<br>&nbsp;
</p>
<p>
This noise generator is based on a function that generates
a random real number uniformely in the semi-open range [0.0, 1.0). 
The function uses the standard Wichmann-Hill generator, 
combining three pure multiplicative congruential generators of 
modulus 30269, 30307 and 30323. Its period (how many numbers it 
generates before repeating the sequence exactly) is 6,953,607,871,644. 
While of much higher quality than the rand() function supplied by 
most C libraries, the theoretical properties are much the same 
as for a single linear congruential generator of large modulus.
</p>
</html>"));

protected
  outer SampleClock sampleClock "Global options";
  Internal.DiscreteNoise discretePart(
    y_min=y_min,
    y_max=y_max,
    firstSeed=firstSeed,
    sampleFactor=sampleFactor) if not continuous "Discrete noise";
equation
  if continuous then
    y = 0.0;
  end if;
connect(y,discretePart.y);

end Noise;
