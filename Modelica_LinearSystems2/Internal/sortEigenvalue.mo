within Modelica_LinearSystems2.Internal;
function sortEigenvalue
  "Sort elements of Eigenvalue-record depending on the eigenvalues"
  extends Modelica.Icons.Function;

  import Modelica_LinearSystems2.Internal.Eigenvalue;
  import Modelica_LinearSystems2.Math.Complex;

  input Eigenvalue ev[:] "Vector to be sorted";
  input Boolean ascending = true
    "= true if ascending order, otherwise descending order";
  input Boolean sortFrequency=true
    "True, if sorting is first for imaginary then for real value, otherwise sorting is for absolute value";
  output Eigenvalue sorted_ev[size(ev,1)] = ev "Sorted vector";
  output Integer indices[size(ev,1)] = 1:size(ev,1) "sorted_ev = ev[indices]";

  /* shellsort algorithm; should be improved later */
protected
  Integer gap;
  Integer i;
  Integer j;
  Eigenvalue wev;
  Integer wi;
  Integer nev = size(ev,1);
  Boolean swap;
  Integer k1;
  Integer k2;
algorithm
  gap := div(nev,2);

  while gap > 0 loop
     i := gap;
     while i < nev loop
        j := i-gap;
        if j>=0 then
           k1 := j+1;
           k2 := j + gap + 1;
           if sortFrequency then
              if ascending then
                 swap := abs(sorted_ev[k1].ev.im) >  abs(sorted_ev[k2].ev.im) or
                         abs(sorted_ev[k1].ev.im) == abs(sorted_ev[k2].ev.im) and
                         (sorted_ev[k1].ev.re  > sorted_ev[k2].ev.re or
                          sorted_ev[k1].ev.re  == sorted_ev[k2].ev.re and sorted_ev[k1].ev.im < sorted_ev[k2].ev.im);
              else
                 swap := abs(sorted_ev[k1].ev.im) <  abs(sorted_ev[k2].ev.im) or
                         abs(sorted_ev[k1].ev.im) == abs(sorted_ev[k2].ev.im) and
                         (sorted_ev[k1].ev.re  < sorted_ev[k2].ev.re or
                          sorted_ev[k1].ev.re  == sorted_ev[k2].ev.re and sorted_ev[k1].ev.im < sorted_ev[k2].ev.im);
              end if;
           else
              if ascending then
                 swap :=Modelica.ComplexMath.'abs'(sorted_ev[k1].ev) > Modelica.ComplexMath.'abs'(sorted_ev[k2].ev);
              else
                 swap :=Modelica.ComplexMath.'abs'(sorted_ev[k1].ev) < Modelica.ComplexMath.'abs'(sorted_ev[k2].ev);
              end if;
           end if;
        else
           swap := false;
        end if;

        while swap loop
           wev := sorted_ev[j+1];
           wi := indices[j+1];
           sorted_ev[j+1] := sorted_ev[j+gap+1];
           sorted_ev[j+gap+1] := wev;
           indices[j+1] := indices[j+gap+1];
           indices[j+gap+1] := wi;
           j := j - gap;
           if j >= 0 then
              k1 := j+1;
              k2 := j + gap + 1;
              if sortFrequency then
                 if ascending then
                    swap := abs(sorted_ev[k1].ev.im) >  abs(sorted_ev[k2].ev.im) or
                            abs(sorted_ev[k1].ev.im) == abs(sorted_ev[k2].ev.im) and
                            (sorted_ev[k1].ev.re  > sorted_ev[k2].ev.re or
                             sorted_ev[k1].ev.re  == sorted_ev[k2].ev.re and sorted_ev[k1].ev.im < sorted_ev[k2].ev.im);
                 else
                    swap := abs(sorted_ev[k1].ev.im) <  abs(sorted_ev[k2].ev.im) or
                            abs(sorted_ev[k1].ev.im) == abs(sorted_ev[k2].ev.im) and
                            (sorted_ev[k1].ev.re  < sorted_ev[k2].ev.re or
                             sorted_ev[k1].ev.re  == sorted_ev[k2].ev.re and sorted_ev[k1].ev.im < sorted_ev[k2].ev.im);
                 end if;
              else
                 if ascending then
                    swap :=Modelica.ComplexMath.'abs'(sorted_ev[k1].ev) > Modelica.ComplexMath.'abs'(sorted_ev[k2].ev);
                 else
                    swap :=Modelica.ComplexMath.'abs'(sorted_ev[k1].ev) < Modelica.ComplexMath.'abs'(sorted_ev[k2].ev);
                 end if;
              end if;
           else
              swap := false;
           end if;
        end while;
        i := i + 1;
     end while;
     gap := div(gap,2);
  end while;
  annotation (Documentation(info="<html>
<h4>Syntax</h4>
<blockquote><pre>
           sorted_v = Vectors.<b>sort</b>(v);
(sorted_v, indices) = Vectors.<b>sort</b>(v, ascending=true);
</pre></blockquote>
<h4>Description</h4>
<p>
Function <b>sort</b>(..) sorts a Real vector v
in ascending order and returns the result in sorted_v.
If the optional argument &quot;ascending&quot; is <b>false</b>, the vector
is sorted in descending order. In the optional second
output argument the indices of the sorted vector with respect
to the original vector are given, such that sorted_v = v[indices].
</p>
<h4>Example</h4>
<blockquote><pre>
  (v2, i2) := Vectors.sort({-1, 8, 3, 6, 2});
       -> v2 = {-1, 2, 3, 6, 8}
          i2 = {1, 5, 3, 4, 2}
</pre></blockquote>
</html>"));
end sortEigenvalue;
