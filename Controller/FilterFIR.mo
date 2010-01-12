within Modelica_LinearSystems2.Controller;
block FilterFIR "Discrete finite impulse response low or high pass filter"
  import Modelica_LinearSystems2.Controller.Types.FIRspec;
  import Modelica_LinearSystems2.Controller.Types.BlockType;
  extends Interfaces.PartialSISO_equality;
  parameter Modelica_LinearSystems2.Controller.Types.FIRspec specType=Modelica_LinearSystems2.Controller.Types.FIRspec.MeanValue
    "Specification type of FIR filter" annotation(Dialog(enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous));
  parameter Integer L(min=2) = 2 "Length of mean value filter" annotation(Dialog(group="Mean value filter",enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.MeanValue));
  parameter Modelica_LinearSystems2.Types.FilterType filterType=
      Modelica_LinearSystems2.Types.FilterType.LowPass "Type of filter" 
                            annotation(Dialog(group="FIR filter design",enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.Window));
  parameter Integer order(min=1) = 2 "Order of filter" annotation(Dialog(group="FIR filter design",enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.Window));
  parameter Modelica.SIunits.Frequency f_cut=1 "Cut-off frequency" annotation(Dialog(group="FIR filter design",enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.Window));
  parameter Types.Window window=Modelica_LinearSystems2.Controller.Types.Window.Rectangle
    "Type of window" annotation(Dialog(group="FIR filter design",enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.Window));
  parameter Real beta=2.12 "Beta-Parameter for Kaiser-window" 
    annotation(Dialog(group="FIR filter design",enable=blockType<>BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.Window and window==Modelica_LinearSystems2.Controller.Types.Window.Kaiser));
  parameter Real a[:]={1,1} "FIR filter coefficients" annotation(Dialog(group="FIR filter defined by coefficient vector",enable=blockType<>Modelica_LinearSystems2.Controller.Types.BlockType.Continuous and specType==Modelica_LinearSystems2.Controller.Types.FIRspec.Coefficients));

  annotation (defaultComponentName="filter",Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics
        ={
        Polygon(
          points={{-82,88},{-90,66},{-74,66},{-82,86},{-82,88}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-82,76},{-82,-92}}, color={192,192,192}),
        Polygon(
          points={{88,-82},{66,-74},{66,-90},{88,-82}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-92,-82},{80,-82}}, color={192,192,192}),
        Text(
          extent={{-14,32},{92,-10}},
          lineColor={192,192,192},
          textString="FIR"),
        Rectangle(
          extent={{-82,-82},{-16,4}},
          lineColor={160,160,164},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Line(points={{-82,30},{-70,30},{-50,28},{-30,20},{-24,16},{-20,12},{-16,
              6},{-12,-4},{-2,-46},{2,-64},{4,-82}}, color={0,0,127}),
        Line(points={{4,-82},{6,-64},{10,-56},{14,-56},{18,-60},{20,-66},{22,-82}}, 
            color={0,0,127}),
        Line(points={{22,-80},{22,-78},{22,-72},{24,-66},{26,-64},{30,-64},{34,
              -66},{36,-70},{38,-78},{38,-82},{38,-74},{40,-68},{42,-66},{46,-66},
              {48,-68},{50,-72},{52,-78},{52,-82},{52,-78},{54,-70},{56,-68},{
              60,-68},{64,-72},{66,-76},{66,-78},{66,-80},{66,-82}}, color={0,0,
              127}),
        Text(
          extent={{-72,86},{98,52}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}));

protected
  parameter Real a2[:]=Internal.FIR_coefficients(
      specType,
      L,
      filterType,
      order,
      f_cut,
      sampleClock.sampleTime*sampleFactor,
      window,
      beta,
      a) if  not continuous;
  Internal.DiscreteFIR discretePart(
    sampleFactor=sampleFactor,
    a=a2) if  not continuous "FIR realization";
equation
 assert(f_cut<=1/(2*sampleClock.sampleTime*sampleFactor),"The cut-off frequency f_cut may not be greater than half the sample frequency (Nyquist frequency), i.e. f_cut <= " + String(1/2/sampleClock.sampleTime*sampleFactor) + " but is "+String(f_cut));
  if continuous then
    y = u;
  end if;
connect(u,discretePart.u);
connect(y,discretePart.y);

end FilterFIR;
