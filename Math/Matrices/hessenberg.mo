within Modelica_LinearSystems2.Math.Matrices;
function hessenberg
  "Compute an upper Hessenberg matrix by repeatedly applicated householder similarity transformation"
  import Modelica;
  import Modelica_LinearSystems2.Math.Matrices;
  import Modelica_LinearSystems2.Math.Vectors;

  input Real A[size(A, 1),size(A, 2)];

  output Real Ht[size(A, 1),size(A, 2)];

algorithm
  Ht := Modelica_LinearSystems2.Math.Matrices.toUpperHessenberg(A, 1, size(A, 1));

  annotation (Documentation(info="<html>
  
   <h4>Syntax</h4>
<blockquote><pre>
         H = Matrices<b>hessenberg</b>(A);
 </pre></blockquote>
<h4>Description</h4>
Function <b>toUpperHessenberg</b> computes the Hessenberg matrix of matrix <b>A</b> by repetitive application of Householder similarity transformation <b>Q</b>'_i*<b>A</b>_i*<b>Q</b>_i
The elementary transformations can be subsumed under
<blockquote> 
<pre>
   <b>H</b> = <b>Q</b>'*<b>A</b>*<b>Q</b>
</pre>
<p>
with
<blockquote> 
<pre>
   <b>Q</b> = <b>Q</b>_1*<b>Q</b>_2*...*<b>Q</b>_n
</pre>
<p>

<h4>Example</h4>
<blockquote><pre>
 A  = [1, 2,  3;
       6, 5,  4;
       1, 0,  0]; 

 H = toUpperHessenberg(A);

  results in:
  
 H = [1.0,  -2.466,  2.630;
     -6.083, 5.514, -3.081;
      0.0,   0.919, -0.514]
      
</pre></blockquote> 

</html>"));
end hessenberg;
