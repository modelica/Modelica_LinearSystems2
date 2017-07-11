within Modelica_LinearSystems2.Controllers.Internal;
block DiscreteFIR "Realization of FIR filter"
  extends Interfaces.PartialDiscreteSISO_equality;
  parameter Real a[:]={1,1} "Coefficients of FIR filter";
protected
  parameter Integer n = size(a, 1) - 1 annotation(HideResult=true);
  discrete Real x[n] annotation(HideResult=true);
  discrete Real sum[n] annotation(HideResult=true);
  discrete Real y_sampled "Sampled output" annotation(HideResult=true);
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
  //u = pre(u);
  x = pre(x);
  //sum = pre(sum);

  annotation (
    Documentation(info="<html>
</html>"), Icon(graphics={Line(
          points={{86,-84},{-94,-84}},
          color={175,175,175},
          smooth=Smooth.None),
        Line(
          points={{-84,76},{-84,-92}},
          color={175,175,175},
          smooth=Smooth.None),
        Polygon(
          points={{-84,90},{-92,68},{-76,68},{-84,90}},
          lineColor={175,175,175},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{90,-84},{68,-92},{68,-76},{90,-84}},
          lineColor={175,175,175},
          smooth=Smooth.None,
          fillColor={175,175,175},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-84,-84},{-18,4}},
          lineColor={175,175,175}),
        Line(
          points={{-84,28},{-72,28},{-52,26},{-32,18},{-26,14},{-22,10},{-18,4},
              {-14,-6},{-4,-48},{0,-66},{2,-84}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{2,-84},{4,-66},{8,-58},{12,-58},{16,-62},{18,-68},{20,-84}},
          color={0,0,127},
          smooth=Smooth.None),
        Line(
          points={{20,-82},{20,-80},{20,-74},{22,-68},{24,-66},{28,-66},{32,-68},
              {34,-72},{36,-80},{36,-84},{36,-76},{38,-70},{40,-68},{44,-68},{46,
              -70},{48,-74},{50,-80},{50,-84},{50,-80},{52,-72},{54,-70},{58,-70},
              {62,-74},{64,-78},{64,-80},{64,-82},{64,-84}},
          color={0,0,127},
          smooth=Smooth.None),
        Text(
          extent={{-20,60},{80,30}},
          lineColor={95,95,95},
          textString="FIR")}));
end DiscreteFIR;
