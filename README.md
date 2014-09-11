# Modelica_LinearSystems2

Free library providing different representations of linear, time invariant differential and difference equation systems, as well as typical operations on these system descriptions.
Additionally, data structures and operations for Complex numbers and for Polynomials are provided. These are utilized above, but are, of course, also useful for other purposes.

## Library description

The `Modelica_LinearSystems2` library is a Modelica package providing different representations of linear, time invariant differential and difference equation systems.

For example, record StateSpace defines a linear time invariant differential equation system in state space form:

        der(x) = A * x + B * u
            y  = C * x + D * u

Operators are overloaded to work conveniently with these system descriptions in an interactive environment, e.g., to multiply transfer functions or to operate on complex numbers.
About 180 functions are provided to operate on these data structures, e.g., to compute eigen values, zeros, step responses, to design pole-placement and LQG controllers, to plot step responses, frequency responses, eigen values, to convert between different description forms, or to generate a linear system description by linearization of a Modelica model.

Furthermore, in sublibrary Controller about 20 input/output blocks of linear systems are provided that are based on the different representation forms, e.g., PID, StateSpace, Filter blocks. A unique feature of these blocks is that it is very convenient to quickly switch between a continuous and a discrete block representation. Also, templates are provide to quickly built-up standard controller structures.

It is planned to include this library in a future version of the Modelica Standard Library.

*Note, the library is not backwards compatible to the previous beta version 0.95, called "Modelica_LinearSystems", which was shipped with previous versions of Dymola. Since the differences are too large, no conversion scripts are provided, but different library names are used.*


## Current release

Download [Modelica_LinearSystems2 v2.3.2 (2014-09-11)](../../archive/v2.3.2.zip)

#### Release notes
*  [Version v2.3.2 (2014-09-11)](../../archive/v2.3.2.zip)
    * This version requires the Modelica Standard Library 3.2.1
    * Improvements performed in version:
        * All Bode diagrams (in all representation forms) can be optionally plotted as magnitude in dB over angular frequency in rad/s, instead of the default to use magnitude over frequency in Hz.
    * Bug fixes performed in version:
        * When computing the gain of the ZerosAndPoles object, a better algorithm is used (the previous one could give bad results if there are large zeros or poles with positive Re-values).
        * Enumeration1/Enumeration2 errors corrected (issued as warning due to stricter checking by Dymola 2015 FD01).
        * Some Enumeration/Integer errors corrected (issued as warning due to stricter checking by Dymola 2015 FD01).
        * Some Plot functions have been called without providing record input arguments. This has been fixed by providing a default record in the function definitions.
*  [Version v2.3.1 (2013-10-01)](../../archive/v2.3.1.zip)
    * This version requires the Modelica Standard Library 3.2.1
        * Backward compatible release to the previous 2.x versions.
*  [Version v2.3 (2012-08-30)](../../archive/v2.3.zip)
    * Backward compatible release to the previous 2.x versions. It contains the following main improvements:
       * New package `Modelica_LinearSystems2.ModelAnalysis` that contains several functions to linearize a model and perform a selected linear analysis operation.
       * New plot functions in package `Modelica_LinearSystems2.Utilities.Plot` to plot parameterized curvces, as well as to plot a root locus of a model, by linearizing a model for a set of selected parameter values.
       * A new function `Modelica_LinearSystems2.Utilities.Import.linearize2` to (a) set parameters of a model, (b) linearize the model and (c) return a StateSpace object (the existing Import.linearize function does only allow to linearize around the default parameter settings and the function returns the A,B,C,D matrices and therefore it is not possible to utilize directly the many functions operating on StateSpace objects).

*  [Version v2.0 (2009-09-02)](../../archive/v2.0.zip)
     * New library based on library Modelica_LinearSystems2 (version 0.95) but is not backwards compatible to this library due to many changes.
     * `Modelica_LinearSystems2` is in principal based on standard Modelica 3.1. However:
       * The plotting functions and the linearization use the Dymola API. This is separated out in package `Modelica_LinearSystems2.Utilities`. The plan is to move this (always tool-dependent package) in to the ModelicaServices package. This is documented under: `Modelica_LinearSystems2.UsersGuide.Requirements`.

## License

This Modelica package is free software and the use is completely at your own risk;
it can be redistributed and/or modified under the terms of the [Modelica License 2](https://modelica.org/licenses/ModelicaLicense2).

## Development and contribution
Release manager: [Martin Otter](http://www.robotic.dlr.de/Martin.Otter)

You may report any issues by using the [Modelica Issue Tracker](https://trac.modelica.org/Modelica/newticket?component=_Modelica_LinearSystems2).

## Acknowledgement
Partial financial support by BMBF (BMBF Förderkennzeichen: 01IS07022) for this work with-in the ITEA2 project [EUROSYSLIB](https://modelica.org/publications/newsletters/2009-1/index_html#eurosyslib) is highly appreciated.
