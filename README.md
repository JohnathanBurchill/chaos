# C CHAOS

This is an implementation of the CHAOS-7 core and static magnetic field models in C for calculating external magnetic field residuals from the [ESA Swarm mission](https://www.esa.int/Applications/Observing_the_Earth/FutureEO/Swarm), primarily for investigations of auroral and ionospheric electrodynamics.

The [CHAOS-7 series of magnetic field models](http://www.spacecenter.dk/files/magnetic-models/CHAOS-7/) are developed and maintained by [DTU Space](https://www.space.dtu.dk/english/research/scientific_data_and_models/magnetic_field_models).

## Swarm 50 Hz residual field

Residual field estimation trades a little bit of accuracy for a lot of speed. Core and crustal (static) magnetic field values are linearly interpolated from control points every 4 s from either the 1 Hz or 50 Hz Swarm MAGx dataset. The model and residual fields are stored in a NASA CDF file.

On a 2022 desktop running GNU/Linux, a daily 50 Hz MAG file takes about 25 s using a single process. This does not include the time it takes to get the unarchived MAGx CDF file onto the local hard drive from the ESA server. Those measurements are available from the ESA Swarm Data Access portal at [1 Hz](https://swarm-diss.eo.esa.int/#swarm%2FLevel1b%2FLatest_baselines%2FMAGx_LR) and [50 Hz](https://swarm-diss.eo.esa.int/#swarm%2FLevel1b%2FLatest_baselines%2FMAGx_HR).

## Magnetic field line of force

The program `tracechaos` traces a magnetic field line of force given a date and an initial position. 

The program `themis_asi_fieldlines` traces field-lines from an assumed auroral emission altitude to a target altitude for the corners of the THEMIS ASI pixels for a given ground station. Calibration files can be generated from THEMIS ASI Level 1 imagery using [AllSkyCameraCal](https://github.com/JohnathanBurchill/AllSkyCameraCal).

## Dependencies

[GNU/Linux](https://www.getgnulinux.org/en/linux), [CMake](https://cmake.org), the [GNU Compiler Collection](https://gcc.gnu.org), the [GNU Scientifc Library](https://www.gnu.org/software/gsl/), and the [NASA CDF library](https://cdf.gsfc.nasa.gov).

## License

Copyright (C) 2022   Johnathan K Burchill

Developed at the University of Calgary (2022) with support of the Canadian Space Agency. 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
