within Modelica_LinearSystems2.Controller.Internal;
block DiscreteFIR "Realization of FIR filter"
  extends Interfaces.PartialDiscreteSISO_equality;
  parameter Real a[:]={1,1} "Coefficients of FIR filter";
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
  annotation (
    Documentation(info="<HTML>
</HTML>
"));
end DiscreteFIR;
