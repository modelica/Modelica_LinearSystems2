within Modelica_LinearSystems2.Controller;
block Filter
  "Continuous or discretized analog low or high pass IIR-filter (CriticalDamping/Bessel/Butterworth/Chebyshev)"
  import Modelica_LinearSystems2.Utilities.Types;
  extends ZerosAndPoles(final system=system3);
  parameter Types.AnalogFilter analogFilter=Modelica_LinearSystems2.Utilities.Types.AnalogFilter.CriticalDamping "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
  parameter Types.FilterType filterType=Modelica_LinearSystems2.Utilities.Types.FilterType.LowPass "Type of filter (LowPass/HighPass)";
  parameter Integer order(min=1) = 2 "Order of filter";
  parameter Modelica.Units.SI.Frequency f_cut=1 "Cut-off frequency";
  parameter Real gain=1.0
    "Gain (= amplitude of frequency response at zero frequency)";
  parameter Boolean normalized=true
    "True, if amplitude at f_cut decreases/increases 3 db (for low/high pass filter), otherwise unmodified filter";
  parameter Real A_ripple(unit="dB") = 0.5
    "Pass band ripple for Chebyshev filter (otherwise not used)" annotation(Dialog(enable=analogFilter == Modelica_LinearSystems2.Utilities.Types.AnalogFilter.Chebyshev));

protected
  parameter Modelica_LinearSystems2.ZerosAndPoles.Internal.ZerosAndPoles
    system2 = Modelica_LinearSystems2.ZerosAndPoles.Internal.filter(
      analogFilter=analogFilter,
      filterType=filterType,
      order=order,
      f_cut=f_cut,
      gain=gain,
      A_ripple=A_ripple,
      normalized=normalized) "Filter" annotation(HideResult=true);

  parameter Modelica_LinearSystems2.ZerosAndPoles system3(
    k=system2.k,
    n1=system2.n1,
    n2=system2.n2,
    d1=system2.d1,
    d2=system2.d2);

  annotation (
    defaultComponentName="filter",
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}), graphics={
        Rectangle(
          extent={{-96,96},{96,-96}},
          lineColor={230,230,255},
          fillColor={230,230,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,80},{-80,-88}}, color={192,192,192}),
        Polygon(
          points={{-80,92},{-88,70},{-72,70},{-80,90},{-80,92}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Line(points={{-90,-78},{82,-78}}, color={192,192,192}),
        Polygon(
          points={{90,-78},{68,-70},{68,-86},{90,-78}},
          lineColor={192,192,192},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-66,60},{88,22}},
          lineColor={192,192,192},
          textString="filter"),
        Text(
          extent={{-136,-104},{164,-134}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="f_cut=%f_cut"),
        Line(points={{22,10},{14,18},{6,22},{-12,28},{-80,28}}, color={0,0,127}),
        Rectangle(
          extent={{-80,-78},{22,10}},
          lineColor={160,160,164},
          fillColor={255,255,255},
          fillPattern=FillPattern.Backward),
        Line(points={{22,10},{30,-2},{36,-20},{40,-32},{44,-58},{46,-78}}),
        Text(
          extent={{-70,96},{98,66}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%sampleFactor")}),
    Documentation(info="<html>
<p>
For details of the filter characteristics, see
<a href=\"modelica://Modelica_LinearSystems2.ZerosAndPoles.Design.filter\">ZerosAndPoles.Design.filter</a>.
</p>
</html>"));

end Filter;
