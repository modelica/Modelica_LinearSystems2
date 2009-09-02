within Modelica_LinearSystems2.Controller.Internal;
block DiscreteSampler "Sample the input signal"
  extends Interfaces.PartialDiscreteSISO_equality;

  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon(
      Ellipse(extent=[-25, -10; -45, 10], style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Ellipse(extent=[45, -10; 25, 10], style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Line(points=[-100, 0; -45, 0], style(color=74, rgbcolor={0,0,127})),
      Line(points=[45, 0; 100, 0], style(color=74, rgbcolor={0,0,127})),
      Line(points=[-35, 0; 30, 35], style(color=74, rgbcolor={0,0,127}))),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram(
      Line(points=[-100, 0; -60, 0]),
      Line(points=[60, 0; 100, 0]),
      Ellipse(extent=[-25, -10; -45, 10], style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Ellipse(extent=[45, -10; 25, 10], style(
          color=74,
          rgbcolor={0,0,127},
          fillColor=7,
          rgbfillColor={255,255,255})),
      Line(points=[-100, 0; -45, 0], style(color=74, rgbcolor={0,0,127})),
      Line(points=[45, 0; 100, 0], style(color=74, rgbcolor={0,0,127})),
      Line(points=[-35, 0; 30, 35], style(color=74, rgbcolor={0,0,127}))),
    Documentation(info="<HTML>
</HTML>
"));
equation
  when {initial(), sampleTrigger} then
      u_sampled = u;
  end when;

  y = u_sampled;
end DiscreteSampler;
