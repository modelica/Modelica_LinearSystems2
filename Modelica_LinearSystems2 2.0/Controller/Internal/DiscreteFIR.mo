within Modelica_LinearSystems2.Controller.Internal;
block DiscreteFIR "Realization of FIR filter"
  extends Interfaces.PartialDiscreteSISO_equality;
  parameter Real a[:]={1,1} "Coefficients of FIR filter";
  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon(
      Line(points=[-84,76; -84,-92],   style(color=8)),
      Polygon(points=[-84,88; -92,66; -76,66; -84,86; -84,88],      style(
          color=8,
          fillColor=8,
          fillPattern=1)),
      Rectangle(extent=[-84,-82; -18,4],   style(
          color=9,
          fillColor=7,
          fillPattern=8)),
      Line(points=[-94,-82; 78,-82],   style(color=8)),
      Polygon(points=[86,-82; 64,-74; 64,-90; 86,-82],     style(
          color=8,
          fillColor=8,
          fillPattern=1)),
      Text(
        extent=[-22,82; 84,40],
        style(color=8),
        string="FIR"),
      Line(points=[-84,30; -72,30; -52,28; -32,20; -26,16; -22,12; -18,6; -14,
              -4; -4,-46; 0,-64; 2,-82],
                                       style(color=74, rgbcolor={0,0,127})),
      Line(points=[2,-82; 4,-64; 8,-56; 12,-56; 16,-60; 18,-66; 20,-82],  style(
            color=74, rgbcolor={0,0,127})),
      Line(points=[20,-80; 20,-78; 20,-72; 22,-66; 24,-64; 28,-64; 32,-66; 34,
              -70; 36,-78; 36,-82; 36,-74; 38,-68; 40,-66; 44,-66; 46,-68; 48,
              -72; 50,-78; 50,-82; 50,-78; 52,-70; 54,-68; 58,-68; 62,-72; 64,
              -76; 64,-78; 64,-80; 64,-82],
                                  style(color=74, rgbcolor={0,0,127}))),
    Window(
      x=0.37,
      y=0.09,
      width=0.52,
      height=0.68),
    Diagram,
    Documentation(info="<HTML>
</HTML>
"));
protected
  parameter Integer n = size(a, 1) - 1 annotation(Hide=true);
  discrete Real x[n] annotation(Hide=true);
  discrete Real sum[n] annotation(Hide=true);
  discrete Real y_sampled "Sampled output" annotation(Hide=true);
equation
  when {initial(), sampleTrigger} then
     u_sampled = u;
     x[1] = pre(u);
     sum[1] = a[2]*x[1];
     x[2:n] = pre(x[1:n - 1]);
     sum[2:n] = a[3:n + 1]*diagonal(x[2:n]) + sum[1:n - 1];
     y_sampled = a[1]*u + sum[n];
  end when;
  y = y_sampled;
initial equation
  u = pre(u);
  x = pre(x);
  sum = pre(sum);
end DiscreteFIR;
