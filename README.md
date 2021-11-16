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

*Note, the library is not backwards compatible to the previous beta version 0.95, called "Modelica_LinearSystems", which was shipped with previous versions of Dymola. Since the differences are too large, no conversion scripts are provided, but different library names are used.*

## Current release

[Modelica_LinearSystems2 v2.4.1+build.2 (2021-11-16)](../../releases/tag/v2.4.1+build.2)

Please note that the library is known to work with Dymola only.

## Older releases

Browse the [Releases](../../releases) page in order to get access to older releases of the `Modelica_LinearSystems2` library.

## License

This Modelica package is free software and the use is completely at your own risk;
it can be redistributed and/or modified under the terms of the [3-Clause BSD License](https://github.com/modelica/ModelicaStandardLibrary/blob/master/LICENSE).

## Development and contribution
The devolopment is organised by [Martin Otter](http://www.robotic.dlr.de/Martin.Otter).

You may report any issues by using the [Issue Tracker](../../issues).

## Acknowledgement
Partial financial support by BMBF (BMBF Förderkennzeichen: 01IS07022) for this work with-in the ITEA2 project [EUROSYSLIB](https://modelica.org/publications/newsletters/2009-1/index_html#eurosyslib) is highly appreciated.
