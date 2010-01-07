within Modelica_LinearSystems2.Controller.Icons;
partial block PartialBlockIcon
  "Basic graphical layout of discrete/continuous block"

protected
                   Boolean cont=true;
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon(
      coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}), 
        graphics={
        Rectangle(
          visible=not cont,
          extent={{-100,100},{100,-100}},
          fillColor={213,255,170},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None),
        Rectangle(
          visible=cont,
          extent={{-100,100},{100,-100}},
          fillColor={230,230,255},
          fillPattern=FillPattern.Solid,
          borderPattern=BorderPattern.Raised,
          pattern=LinePattern.None,
          lineColor={0,0,0}),
        Text(
          extent={{-150,110},{150,150}},
          lineColor={0,0,255},
          textString="%name")}),
    Window(
      x=0.33,
      y=0.33,
      width=0.59,
      height=0.43));

end PartialBlockIcon;
