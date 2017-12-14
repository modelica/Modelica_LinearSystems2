within Modelica_LinearSystems2.WorkInProgress.Tests;
package Filters
  function plotFilter
    import ZP = Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Utilities.Types;
    import Modelica.Utilities.Streams.print;
    import Modelica.Constants.pi;
    import Complex;

    input Modelica_LinearSystems2.Utilities.Types.AnalogFilter analogFilter "Analog filter characteristics (CriticalDamping/Bessel/Butterworth/Chebyshev)";
    input Modelica_LinearSystems2.Utilities.Types.FilterType filterType=Utilities.Types.FilterType.LowPass "Type of filter (LowPass/HighPass)";
     input Modelica.SIunits.Frequency f_cut=1/(2*Modelica.Constants.pi)
      "Cut-off frequency";
     input Real A_ripple(unit="dB") = 0.5
      "Pass band ripple for Chebyshev filter (otherwise not used)";
     input Modelica.SIunits.Frequency f_min=0
      "Band of pass band filter is f_min (-3db) .. f_cut (-3db)";
     input Boolean normalized=true
      "True, if amplitude at f_cut decreases/increases 3 db (for low/high pass filter), otherwise unmodified filter";
  protected
    function getAmplitude "Compute amplitude at f_cut and return it as string"
      import Modelica_LinearSystems2;
      input ZP zp;
      input Modelica.SIunits.Frequency f_cut;
      output String str;
    protected
      Complex c;
      Real A;
      constant Real machEps = 100*Modelica.Constants.eps;
    algorithm
      c :=ZP.Analysis.evaluate(zp, Complex(0, 2*pi*f_cut));
      A :=Modelica.ComplexMath.'abs'(c);
      str :="amplitude(f=" + String(f_cut) + ") = ";
      if Modelica_LinearSystems2.Math.isEqual(A, 10^(-3/20), machEps) then
         str := str + "-3db";
      elseif Modelica_LinearSystems2.Math.isEqual(A, sqrt(2)/2, machEps) then
         str := str + "sqrt(2)/2";
      else
         str := str + String(A);
      end if;
    end getAmplitude;

     function getFilterName "Return the filter name as string"
      input Modelica_LinearSystems2.Utilities.Types.AnalogFilter analogFilter;
        output String str;
     algorithm
      if analogFilter == Modelica_LinearSystems2.Utilities.Types.AnalogFilter.CriticalDamping then
           str :="CriticalDamping";
      elseif analogFilter == Modelica_LinearSystems2.Utilities.Types.AnalogFilter.Bessel then
           str :="Bessel";
      elseif analogFilter == Modelica_LinearSystems2.Utilities.Types.AnalogFilter.Butterworth then
           str :="Butterworth";
      elseif analogFilter == Modelica_LinearSystems2.Utilities.Types.AnalogFilter.Chebyshev then
           str :="Chebyshev";
        else
           str :="unknown";
        end if;
     end getFilterName;

     Real f_0;

     ZP filter1 = ZP.Design.filter(analogFilter, filterType, order=1, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized, f_min = f_min);
     ZP filter2 = ZP.Design.filter(analogFilter, filterType, order=2, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized, f_min = f_min);
     ZP filter3 = ZP.Design.filter(analogFilter, filterType, order=3, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized, f_min = f_min);
     ZP filter4 = ZP.Design.filter(analogFilter, filterType, order=4, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized, f_min = f_min);
     ZP filter5 = ZP.Design.filter(analogFilter, filterType, order=5, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized, f_min = f_min);
     ZP filter6 = ZP.Design.filter(analogFilter, filterType, order=6, f_cut=f_cut, A_ripple=A_ripple, normalized=normalized, f_min = f_min);
  algorithm
     ZP.Plot.bode(filter1);
     ZP.Plot.bode(filter2);
     ZP.Plot.bode(filter3);
     ZP.Plot.bode(filter4);
     ZP.Plot.bode(filter5);
     ZP.Plot.bode(filter6);

     // Check filter properties
     if filterType == Types.FilterType.LowPass then
        print("\nLow pass filter: analogFilter = " + getFilterName(analogFilter));
        print("filter1: dcGain = " + String(ZP.Analysis.dcGain(filter1)) + ", " + getAmplitude(filter1,f_cut));
        print("filter2: dcGain = " + String(ZP.Analysis.dcGain(filter2)) + ", " + getAmplitude(filter2,f_cut));
        print("filter3: dcGain = " + String(ZP.Analysis.dcGain(filter3)) + ", " + getAmplitude(filter3,f_cut));
        print("filter4: dcGain = " + String(ZP.Analysis.dcGain(filter4)) + ", " + getAmplitude(filter4,f_cut));
        print("filter5: dcGain = " + String(ZP.Analysis.dcGain(filter5)) + ", " + getAmplitude(filter5,f_cut));
        print("filter6: dcGain = " + String(ZP.Analysis.dcGain(filter6)) + ", " + getAmplitude(filter6,f_cut));

     elseif filterType == Types.FilterType.HighPass then
        print("\nHigh pass filter: analogFilter = " + getFilterName(analogFilter));
        print("filter1: k = " + String(filter1.k) + ", " + getAmplitude(filter1,f_cut));
        print("filter2: k = " + String(filter2.k) + ", " + getAmplitude(filter2,f_cut));
        print("filter3: k = " + String(filter3.k) + ", " + getAmplitude(filter3,f_cut));
        print("filter4: k = " + String(filter4.k) + ", " + getAmplitude(filter4,f_cut));
        print("filter5: k = " + String(filter5.k) + ", " + getAmplitude(filter5,f_cut));
        print("filter6: k = " + String(filter6.k) + ", " + getAmplitude(filter6,f_cut));

     elseif filterType == Types.FilterType.BandPass then
        print("\nBand pass filter: analogFilter = " + getFilterName(analogFilter));
        print("filter1: " + getAmplitude(filter1, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter1, f_min) + ", " +
                            getAmplitude(filter1, f_cut));
        print("filter2: " + getAmplitude(filter2, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter2, f_min) + ", " +
                            getAmplitude(filter2, f_cut));
        print("filter3: " + getAmplitude(filter3, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter3, f_min) + ", " +
                            getAmplitude(filter3, f_cut));
        print("filter4: " + getAmplitude(filter4, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter4, f_min) + ", " +
                            getAmplitude(filter4, f_cut));
        print("filter5: " + getAmplitude(filter5, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter5, f_min) + ", " +
                            getAmplitude(filter5, f_cut));
        print("filter6: " + getAmplitude(filter6, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter6, f_min) + ", " +
                            getAmplitude(filter6, f_cut));

     elseif filterType == Types.FilterType.BandStop then
        print("\nBand stop filter: analogFilter = " + getFilterName(analogFilter));
        print("filter1: " + "dcGain = " + String(ZP.Analysis.dcGain(filter1)) + ", " +
                            getAmplitude(filter1, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter1, f_min) + ", " +
                            getAmplitude(filter1, f_cut));
        print("filter2: " + "dcGain = " + String(ZP.Analysis.dcGain(filter2)) + ", " +
                            getAmplitude(filter2, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter2, f_min) + ", " +
                            getAmplitude(filter2, f_cut));
        print("filter3: " + "dcGain = " + String(ZP.Analysis.dcGain(filter3)) + ", " +
                            getAmplitude(filter3, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter3, f_min) + ", " +
                            getAmplitude(filter3, f_cut));
        print("filter4: " + "dcGain = " + String(ZP.Analysis.dcGain(filter4)) + ", " +
                            getAmplitude(filter4, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter4, f_min) + ", " +
                            getAmplitude(filter4, f_cut));
        print("filter5: " + "dcGain = " + String(ZP.Analysis.dcGain(filter5)) + ", " +
                            getAmplitude(filter5, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter5, f_min) + ", " +
                            getAmplitude(filter5, f_cut));
        print("filter6: " + "dcGain = " + String(ZP.Analysis.dcGain(filter6)) + ", " +
                            getAmplitude(filter6, sqrt(f_min*f_cut)) + ", " +
                            getAmplitude(filter6, f_min) + ", " +
                            getAmplitude(filter6, f_cut));
     end if;

   annotation(__Dymola_interactive=true);
  end plotFilter;

  function plotFilter2
    import ZP = Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Utilities.Types;
    input Modelica_LinearSystems2.Utilities.Types.FilterType filterType "Type of filter (LowPass/HighPass)";
    input Modelica.SIunits.Frequency f_cut=3 "Cut-off frequency";
    input Real A_ripple(unit="dB") = 0.5
      "Pass band ripple for Chebyshev filter (otherwise not used)";
    input Modelica.SIunits.Frequency f_min=0
      "Band of pass band filter is f_min (-3db) .. f_cut (-3db)";
    input Boolean normalized=true
      "True, if amplitude at f_cut decreases/increases 3 db (for low/high pass filter), otherwise unmodified filter";
  algorithm
    plotFilter(Types.AnalogFilter.CriticalDamping, filterType, f_cut, A_ripple, f_min, normalized);
    plotFilter(Types.AnalogFilter.Bessel,          filterType, f_cut, A_ripple, f_min, normalized);
    plotFilter(Types.AnalogFilter.Butterworth,     filterType, f_cut, A_ripple, f_min, normalized);
    plotFilter(Types.AnalogFilter.Chebyshev,       filterType, f_cut, A_ripple, f_min, normalized);

    annotation(__Dymola_interactive=true);
  end plotFilter2;

  function plotFilter3
    import ZP = Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica_LinearSystems2.Utilities.Types;
    input Modelica.SIunits.Frequency f_cut=3 "Cut-off frequency";
    input Real A_ripple(unit="dB") = 0.5
      "Pass band ripple for Chebyshev filter (otherwise not used)";
    input Modelica.SIunits.Frequency f_min=2
      "Band of pass band filter is f_min (-3db) .. f_cut (-3db)";
    input Boolean normalized=true
      "True, if amplitude at f_cut decreases/increases 3 db (for low/high pass filter), otherwise unmodified filter";
  algorithm
    //plotFilter2(Types.FilterType.LowPass,  f_cut, A_ripple, f_min, normalized);
    //plotFilter2(Types.FilterType.HighPass, f_cut, A_ripple, f_min, normalized);
    // plotFilter2(Types.FilterType.BandPass, f_cut, A_ripple, f_min, normalized);
    plotFilter2(Types.FilterType.BandStop, f_cut, A_ripple, f_min, normalized);

    annotation(__Dymola_interactive=true);
  end plotFilter3;

  function compareBaseFiltersWithTietzeSchenk
    "Compare normalized base filters with the table of Tietze/Schenk Halbleiterschaltungstechnik"
    import Modelica_LinearSystems2.Utilities.Types;
    import ZP = Modelica_LinearSystems2.ZerosAndPoles;
    import Modelica.Utilities.Streams.print;
    import Complex;
    input String outputFile = "";
  protected
    constant Real machEps = 100*Modelica.Constants.eps;
    constant Real eps = 0.001;
    ZP zp;
    Integer maxOrder = 5;
    Boolean evenOrder "True, if even filter order, otherwise uneven";
    Real c;
    Real k;
    Boolean gainIsOne;

    function printCoefficients
      "Transform coefficients to Tietze/Schenk and print them"
      /*  Tietze/Schenk:  (a1*p + 1); (b2*p^2 + a2*p + 1)
        ZerosAndPoles:  (p + a3)  ; (p^2 + b3*p + a3)
        Therefore:      a1 = 1/a3 ; b2 = 1/a3, a2 = b3/a3
    */
      input ZP zp;
      input String outputFile = "";
    protected
      Real k;
      Boolean gainIsOne;
    algorithm
      k :=ZP.Analysis.dcGain(zp);
      gainIsOne := Modelica_LinearSystems2.Math.isEqual(k, 1.0, machEps);
      if not gainIsOne then
         print("!!! Gain of base filter is wrong (should be one, but is " + String(k) + ")", outputFile);
      end if;

      for i in 1:size(zp.d1,1) loop
         print("  a = " + String(1/zp.d1[i],format="6.5f") + ", b = 0", outputFile);
      end for;

      for i in 1:size(zp.d2,1) loop
         print("  a = " + String(zp.d2[i,1]/zp.d2[i,2],format="6.5f") +
               ", b = " + String(1/zp.d2[i,2],format="6.5f"), outputFile);
      end for;
    end printCoefficients;

    function getAmplitude "Compute amplitude at w=1 and return it as string"
      import Modelica_LinearSystems2;
      input ZP zp;
      output String str;
    protected
      Complex c;
      Real A;
    algorithm
      c :=ZP.Analysis.evaluate(zp, Complex(0, 1.0));
      A :=Modelica.ComplexMath.'abs'(c);

      if Modelica_LinearSystems2.Math.isEqual(A, 10^(-3/20), machEps) then
         str :="amplitude(w=1) = -3db";
      elseif Modelica_LinearSystems2.Math.isEqual(A, sqrt(2)/2, machEps) then
         str :="amplitude(w=1) = sqrt(2)/2";
      else
         str :="amplitude(w=1) = " + String(A);
      end if;
    end getAmplitude;
  algorithm
    print("\n" +
          "... The following values should be identical to the tables in Abb. 13.14\n"+
          "... of Tietze/Schenk 2002, pp. 828-834", outputFile);

    // Critical damping
    print("\nCriticalDamping filter (all even coefficients are identical):", outputFile);

    for i in 1:maxOrder loop
       zp :=Modelica_LinearSystems2.ZerosAndPoles.Internal.baseFilter(
                                 Types.AnalogFilter.CriticalDamping, order=i);
       print("\n  order = " + String(i) + ", " + getAmplitude(zp), outputFile);

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
       assert(gainIsOne, "CriticalDamping base filter is wrong (3)");

       /* Check coefficients of first and second order transfer functions
        Tietze/Schenk:  (a1*p + 1); (b2*p^2 + a2*p + 1)
        ZerosAndPoles:  (p + a3)  ; (p + a3)^2 = p^2 + 2*a3*p + a3^2)
        Therefore:       a1 = 1/a3; a2 = 2/a3, b2 = 1/a3^2
     */
       evenOrder :=mod(i, 2) == 0;
       if i==1 then
          print("  a = " + String(1/zp.d1[1]), outputFile);

       elseif evenOrder then
          print("  a = " + String(2/zp.d1[1]) + ", b = " + String(1/zp.d1[1]^2), outputFile);

       else
          print("  a = " + String(1/zp.d1[1]) + "\n" +
                "  a = " + String(2/zp.d1[1]) + ", b = " + String(1/zp.d1[1]^2), outputFile);

       end if;
    end for;

    // Bessel filter
    print("\nBessel filter:", outputFile);
    for i in 1:maxOrder loop
       zp :=Modelica_LinearSystems2.ZerosAndPoles.Internal.baseFilter(
                                 Types.AnalogFilter.Bessel, order=i);
       print("\n  order = " + String(i)  + ", " + getAmplitude(zp), outputFile);
       printCoefficients(zp, outputFile);
    end for;

    // Butterworth filter
    print("\nButterworth filter:", outputFile);
    for i in 1:maxOrder loop
       zp :=Modelica_LinearSystems2.ZerosAndPoles.Internal.baseFilter(
                                 Types.AnalogFilter.Butterworth, order=i);
       print("\n  order = " + String(i)  + ", " + getAmplitude(zp), outputFile);
       printCoefficients(zp, outputFile);
    end for;

    // Chebyshev filter
    print("\nChebyshev filter (A_ripple = 0.5 db):", outputFile);
    for i in 1:maxOrder loop
       zp :=Modelica_LinearSystems2.ZerosAndPoles.Internal.baseFilter(
                                 Types.AnalogFilter.Chebyshev, A_ripple=0.5, order=i);
       print("\n  order = " + String(i)  + ", " + getAmplitude(zp), outputFile);
       printCoefficients(zp, outputFile);
    end for;

    print("\nChebyshev filter (A_ripple = 3 db):", outputFile);
    for i in 1:maxOrder loop
       zp :=Modelica_LinearSystems2.ZerosAndPoles.Internal.baseFilter(
                                 Types.AnalogFilter.Chebyshev, A_ripple=3, order=i);
       print("\n  order = " + String(i)  + ", " + getAmplitude(zp), outputFile);
       printCoefficients(zp, outputFile);
    end for;
    annotation (Documentation(info="<html>
<p>
This function compares the filters with the ones from
<a href=\"modelica://Modelica_LinearSystems2.UsersGuide.Literature\">[Tietze2002]</a>, pp. 815-852.
</p>

<p>
In tables Abb. 13.14 on pages 828-834 of (Tietze/Schenk 2002), the filter coefficients are
given with respect to normalized filters for a cut-off angular frequency
of 1 rad/s. The normalization is performed in such a way that at the cut-off
frequency the transfer function has an amplitude of -3db (= 10^(-3/20) = 0.7079457..).
In the tables, not the exact -3db value is used but the approximation
sqrt(2)/2 (= 0.707106...). Due to &quot;historical&quot; reasons, function baseFilter
from the Modelica_LinearSystems library uses -3db for Bessel and Chebyshev filters
and sqrt(2)/2 for CriticalDamping and Butterworth filters. Furthermore, the table
gives the values only up to 4 significant digits. For these reasons, in this test
function the comparison is performed up to 3 significant digits.
</p>

</html>"));
  end compareBaseFiltersWithTietzeSchenk;

  function plotAlphaForPassBand
     input Real a = 0.1;
     input Real b = 0.2;
     input Real w = 1;
  protected
     Real z[:]=  0:0.01:2;
     Real res[size(z,1)];
  algorithm
    for i in 1:size(z,1) loop
      res[i] :=z[i]^2 + a^2*w^2*z[i]^2/(1 + z[i])^2 - (2 + b*w^2)*z[i] + 1;
    end for;

    Modelica_LinearSystems2.Utilities.Plot.diagram(
      Modelica_LinearSystems2.Utilities.Plot.Records.Diagram(
          curve={Modelica_LinearSystems2.Utilities.Plot.Records.Curve(
            x=z,
            y=res,
            legend="residue")}));

    annotation(__Dymola_interactive=true);
  end plotAlphaForPassBand;
end Filters;
