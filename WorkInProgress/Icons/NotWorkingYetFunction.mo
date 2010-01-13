within Modelica_LinearSystems2.Icons;
partial function NotWorkingYetFunction
  "Icon for a function which is definitely not working (yet)"

  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,
            -100},{100,100}}), graphics={
        Text(
          extent={{-140,162},{136,102}},
          textString="%name",
          lineColor={0,0,255}),
        Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={255,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-100,100},{100,-100}},
          color={255,0,0},
          smooth=Smooth.None,
          thickness=0.5),
        Line(
          points={{-100,-100},{100,100}},
          color={255,0,0},
          thickness=0.5,
          smooth=Smooth.None)}), Documentation(info="<html>
<p>
This icon is designed for a <b>function</b>
</p>
</html>"));

end NotWorkingYetFunction;
