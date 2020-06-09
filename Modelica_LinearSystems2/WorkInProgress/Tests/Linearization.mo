within Modelica_LinearSystems2.WorkInProgress.Tests;
package Linearization
  model FirstOrder
    extends Modelica.Blocks.Interfaces.SISO;
    Modelica.Blocks.Continuous.FirstOrder firstOrder(k=k, T=T)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}})));
    parameter Real k=1 "Gain";
    parameter Modelica.Units.SI.Time T=1 "Time Constant";
  equation
    connect(firstOrder.u, u) annotation (Line(
        points={{-12,0},{-120,0}},
        color={0,0,127},
        smooth=Smooth.None));
    connect(firstOrder.y, y) annotation (Line(
        points={{11,0},{110,0}},
        color={0,0,127},
        smooth=Smooth.None));
    annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
              -100},{100,100}}), graphics));
  end FirstOrder;
end Linearization;
