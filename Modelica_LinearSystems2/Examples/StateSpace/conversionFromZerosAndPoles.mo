within Modelica_LinearSystems2.Examples.StateSpace;
function conversionFromZerosAndPoles
  "Transform a zeros-and-poles representation into a StateSpace description"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.StateSpace;
  import Modelica_LinearSystems2.ZerosAndPoles;
  import Complex;
  import Modelica.ComplexMath.j;

protected
  ZerosAndPoles zp= ZerosAndPoles({2+0*j}, {1+0*j,2+3*j,2-3*j}, 4);
  StateSpace ss=StateSpace(zp);
public
  output Boolean ok "Standard output";

algorithm
  Modelica.Utilities.Streams.print("zp = " + String(zp));
  Modelica.Utilities.Streams.print("ss = " + String(ss));
  ok := true;

  annotation (
    Documentation(info="<html>
<p>
This example demonstrates the conversion of a zeros-and-poles system into a state space system.
Running this function the following output will be printed containing the input zeroes and poles 
description&nbsp;<code>zp</code> and the resulting output state space description&nbsp;<code>ss</code>.
</p>
<blockquote><pre>
zp = 4*(p - 2) /  ( (p - 1)*(p^2 - 4*p + 13) )
ss = 
  ss.A = 
             x1     x2     x3 
       x1    0      1      0                      
       x2   -13     4      0                      
       x3   -1      0.5    1                      
   
  ss.B = 
             u1 
       x1    0                      
       x2    13                     
       x3    0                      
   
  ss.C = 
             x1     x2     x3 
       y1    0      0      0.615384615385         
   
  ss.D = 
             u1 
       y1    0                      
</pre></blockquote>
</html>"));
end conversionFromZerosAndPoles;
