within Modelica_LinearSystems2.Controllers;
block TransferFunction
  "Continuous or discrete, single input single output transfer function"
  import Modelica_LinearSystems2.Controllers.Types.Init;

  extends Modelica_LinearSystems2.Controllers.Interfaces.PartialSISO2(
      discretePart(
      x_start=x_start,
      y_start={y_start},
      ABCD=Modelica_LinearSystems2.TransferFunction.Conversion.toMatrices(system)));

  parameter Modelica_LinearSystems2.TransferFunction system "Transfer function";
  parameter Real x_start[nx]=zeros(nx) "Initial or guess values of states"
    annotation(Dialog(tab="Advanced options"));
  parameter Real y_start=0
    "Initial value of output (derivatives of y are zero upto nx-1-th derivative)" annotation(Dialog(tab="Advanced options"));
  final parameter Integer nx=size(system.d, 1) - 1 "Number of states x";
  Modelica.Blocks.Interfaces.RealOutput x[nx]
    "State of continuous transfer function";

protected
  parameter Boolean withDelay=false;
  parameter Integer na=size(system.d, 1) annotation(HideResult=true);
  parameter Integer nb=size(system.n, 1) annotation(HideResult=true);
  Real a[na]=system.d "Reverse element order of system.denominator" annotation(HideResult=true);
  Real b[nb]=system.n annotation(HideResult=true);
  Real bb[:]=vector([zeros(max(0, na - nb), 1); b]);
  Real d=bb[1]/a[1];
  Real a_end=if a[end] > 100*Modelica.Constants.eps*sqrt(a*a) then a[end] else 1.0;
  Real x_scaled[size(x, 1)];//=x*a_end "Scaled vector x";

equation
  x_scaled =x*a_end;

  if continuous then
    if nx == 0 then
      y = d*u;
    else
      der(x_scaled[1]) = (-a[2:na]*x_scaled + a_end*u)/a[1];
      der(x_scaled[2:nx]) = x_scaled[1:nx - 1];
      y = ((bb[2:na] - d*a[2:na])*x_scaled)/a_end + d*u;
    end if;
  end if;
  connect(y, discretePart.y[1]);
  connect(x, discretePart.x);
initial equation
  if continuous and nx>0 then
    if init == Init.SteadyState then
      der(x_scaled) = zeros(nx);
    elseif init == Init.InitialState then
      x_scaled = x_start*a_end;
    elseif init == Init.InitialOutput then
      y = y_start;
//      der(x_scaled[1:nx-1]) = zeros(nx -1);
      der(x_scaled[1:nx-(if nx>1 then 2 else 1)]) = zeros(nx -( if nx>1 then 2 else 1));
    end if;
  end if;

  annotation (
    defaultComponentName="transferFunction",
    Icon(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}),
      graphics={
        Text(
          extent={{-90,10},{90,90}},
          textColor={0,0,127},
          textString="n(s)"),
        Line(points={{-80,0},{80,0}}, color={0,0,127}),
        Text(
          extent={{-90,-10},{90,-90}},
          textColor={0,0,127},
          textString="d(s)"),
        Text(
          extent={{-96,-106},{100,-140}},
          textColor={0,0,0},
          textString="%sampleFactor")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        grid={2,2}),
      graphics={
        Line(points={{-100,0},{-60,0}}, color={0,0,255}),
        Rectangle(extent={{-60,60},{60,-60}}, lineColor={0,0,255}),
        Line(points={{40,0},{-40,0}}, color={0,0,0}),
        Text(
          extent={{-55,55},{55,5}},
          textColor={0,0,0},
          textString="n(s)"),
        Text(
          extent={{-55,-5},{55,-55}},
          textColor={0,0,0},
          textString="d(s)"),
        Line(points={{60,0},{100,0}}, color={0,0,255})}),
    Documentation(info="<html>
</html>"));
end TransferFunction;
