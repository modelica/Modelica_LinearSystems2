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
end Filters;
