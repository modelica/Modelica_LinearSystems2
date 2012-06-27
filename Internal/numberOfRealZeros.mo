within Modelica_LinearSystems2.Internal;
function numberOfRealZeros
  "Determine the number of elements of a Complex vector where the imaginary part is zero"
  import Modelica_LinearSystems2.Math.Complex;

  input Complex complexVector[:] "Complex vector";
  output Integer result "Number of elements of v with v.im = 0";

algorithm
  result := 0;
  for i in 1:size(complexVector, 1) loop
    if complexVector[i].im == 0 then
      result := result + 1;
    end if;
  end for;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
     Matrices.<b>numberOfRealZeros</b>(complexVector);
</pre></blockquote>
<h4>Description</h4>
<p>
Function <b>numberOfRealZeros</b>(..) determines the number of
elements of vector \"complexVector\" with vanishing imaginary part,
i.e., complexVector[i].im = 0.
</p>
<h4>Example</h4>
<blockquote><pre>
  // c = {0; 1+2*j; 1-2*j; 2; -3; -1-j; -1+j};
  result = Matrices.numberOfRealZeros(c);
           -> result = 3;
</pre></blockquote>
</html>"));
end numberOfRealZeros;
