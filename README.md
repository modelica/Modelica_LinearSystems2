# Modelica_LinearSystems2

Free library providing different representations of linear, time invariant differential and difference equation systems, as well as typical operations on these system descriptions.

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

Download [Modelica_LinearSystems2 v2.0 (2009-09-02)](../../archive/v2.0.zip)

#### Release notes
*  [Version v2.0 (2009-09-02)](../../archive/v2.0.zip)
 * New library based on library Modelica_LinearSystems2 (version 0.95) but is not backwards compatible to this library due to many changes.
 * `Modelica_LinearSystems2` is in principal based on standard Modelica 3.1. However:
   * The plotting functions and the linearization use the Dymola API. This is separated out in package `Modelica_LinearSystems2.Utilities`. The plan is to move this (always tool-dependent package) in to the ModelicaServices package. This is documented under: `Modelica_LinearSystems2.UsersGuide.Requirements`.

## License

This Modelica package is free software and the use is completely at your own risk;
it can be redistributed and/or modified under the terms of the [Modelica License 2](https://modelica.org/licenses/ModelicaLicense2).

## Development and contribution
Release manager: [Martin Otter](mailto:Martin.Otter@dlr.de)

You may report any issues with using the [Modelica Issue Tracker](https://trac.modelica.org/Modelica/newticket?component=_Modelica_LinearSystems2).

## Acknowledgement
Partial financial support by BMBF (BMBF Förderkennzeichen: 01IS07022) for this work with-in the ITEA2 project EUROSYSLIB is highly appreciated.
