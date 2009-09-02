within Modelica_LinearSystems2.Controller.Interfaces;
partial block PartialBlockIcon
  "Basic graphical layout of discrete/continuous block"

  annotation (
    Coordsys(
      extent=[-100, -100; 100, 100],
      grid=[2, 2],
      component=[20, 20]),
    Icon(Rectangle(extent=[-100, -100; 100, 100], style(
            color=74,
            rgbcolor={0,0,127},
            fillColor=30,
            rgbfillColor={230,230,255},
            fillPattern=11)),
         Text(extent=[-150, 150; 150, 110], string="%name")),
    Window(
      x=0.33,
      y=0.33,
      width=0.59,
      height=0.43));
equation

end PartialBlockIcon;
