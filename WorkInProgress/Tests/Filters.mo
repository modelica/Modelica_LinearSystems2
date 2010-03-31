within Modelica_LinearSystems2.WorkInProgress.Tests;
package Filters
  function PlotFilter
    import ZP = Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Types;
     input Modelica_LinearSystems2.Types.AnalogFilter analogFilter=Types.AnalogFilter.CriticalDamping
      "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
     input Modelica_LinearSystems2.Types.FilterType filterType=Types.FilterType.LowPass
      "Type of filter (LowPass/HighPass)";
     input Modelica.SIunits.Frequency f_cut=3 "Cut-off frequency";
     input Real A_ripple(unit="dB") = 0.5
      "Pass band ripple for Chebyshev filter (otherwise not used)";
     input Boolean normalized=true
      "= true, if amplitude at f_cut decreases/increases 3 db (for low/high pass filter), otherwise unmodified filter";
  protected
     ZP filter1 = ZP.Design.filter(analogFilter, filterType, order=1, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter2 = ZP.Design.filter(analogFilter, filterType, order=2, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter3 = ZP.Design.filter(analogFilter, filterType, order=3, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter4 = ZP.Design.filter(analogFilter, filterType, order=4, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter5 = ZP.Design.filter(analogFilter, filterType, order=5, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter6 = ZP.Design.filter(analogFilter, filterType, order=6, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter7 = ZP.Design.filter(analogFilter, filterType, order=7, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
     ZP filter8 = ZP.Design.filter(analogFilter, filterType, order=8, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized);
  algorithm
     ZP.Plot.bode(filter1, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter2, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter3, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter4, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter5, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter6, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter7, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
     ZP.Plot.bode(filter8, autoRange=false, f_min=0.2, f_max=f_cut*1.2);
  end PlotFilter;

  function CompareBaseFiltersWithTietzeSchenk
    "Compare normalized base filters with the table of Tietze/Schenk Halbleiterschaltungstechnik"
    import Modelica_LinearSystems2.Types;
    import ZP = Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica.Utilities.Streams.print;
    input String outputFile = "";
  protected
    constant Real machEps = 100*Modelica.Constants.eps;
    constant Real eps = 0.001;
    ZP zp;
    Integer maxOrder = 5;
    Boolean evenOrder "= true, if even filter order (otherwise uneven)";
    Real c;
    Real k;
    Boolean gainIsOne;
  algorithm
    print("\n" +
          "... The following values should be identical to the tables in Abb. 13.14\n"+
          "... of Tietze/Schenk 2002, pp. 828-834\n", outputFile);

    // Critical damping
    print("CriticalDamping filter:", outputFile);

    for i in 1:maxOrder loop
       zp :=ZP.Design.baseFilter(Types.AnalogFilter.CriticalDamping, order=i);

       // Check that all coefficients are identical
       assert(size(zp.n1,1) == 0 and
              size(zp.n2,1) == 0 and
              size(zp.d2,1) == 0, "CriticalDamping base filter is wrong (1)");
       c :=zp.d1[1];
       for j in 2:i loop
          assert(zp.d1[i] == c, "CriticalDamping base filter is wrong (2)");
       end for;

       // Check that dc gain is one
       k :=ZP.Analysis.dcGain(zp);
       gainIsOne := Modelica_LinearSystems2.Math.isEqual(k,1.0,machEps);
       assert(gainIsOne,
               "CriticalDamping base filter is wrong (3)");

       /* Check coefficients of first and second order transfer functions
        Tietze/Schenk:  (a1*p + 1); (b2*p^2 + a2*p + 1)
        ZerosAndPoles:  (p + a3)  ; (p + a3)^2 = p^2 + 2*a3*p + a3^2)
        Therefore:       a1 = 1/a3; a2 = 2/a3, b2 = 1/a3^2 
     */
       evenOrder :=mod(i, 2) == 0;
       if i==1 then
          print("order = " + String(i) + ", a1 = " + String(1/zp.d1[1]), outputFile);

       elseif evenOrder then
          print("order = " + String(i) +
                ", a2 = " + String(2/zp.d1[1]) + ", b2 = " + String(1/zp.d1[1]^2), outputFile);

       else
          print("order = " + String(i) + ", a1 = " + String(1/zp.d1[1]) + "\n" +
                "           a2 = " + String(2/zp.d1[1]) + ", b2 = " + String(1/zp.d1[1]^2), outputFile);

       end if;
    end for;

    annotation (Documentation(info="<html>
<p>
This function compares the filters with the ones from
</p>

<dl>
<dt>Tietze U., and Schenk C. (2002):</dt>
<dd> <b>Halbleiter-Schaltungstechnik</b>.
     Springer Verlag, 12. Auflage, pp. 815-852.</dd>
</dl>

<p>
In tables Abb. 13.14 on pages 828-834 of (Tietze/Schenk 2002), the filter coefficients are
given with respect to normalized filters for a cut-off angular frequency
of 1 rad/s. The normalization is performed in such a way that at the cut-off
frequency the transfer function has an amplitude of -3db (= 10^(-3/20) = 0.7079457..).
In the tables, not the exact -3db value is used but the approximation 
sqrt(2)/2 (= 0.707106...). Due to \"historical\" reasons, function baseFilter
from the Modelica_LinearSystems library uses -3db for Bessel and Chebyshev filters
and sqrt(2)/2 for CriticalDamping and Butterworth filters. Furthermore, the table 
gives the values only up to 4 significant digits. For these reasons, in this test
function the comparison is performed up to 3 significant digits.
</p>

</html>"));
  end CompareBaseFiltersWithTietzeSchenk;
end Filters;
