# CHAOS C implementation

This is an implementation of the CHAOS magnetic field model in C. Core and crustal model field values are interpolated every 4 s at times from either the 1 Hz or 50 Hz Swarm MAGx dataset, and residual field vectors are calculated. Results are stored in a NASA CDF file.

On a contemporary desktop machine running GNU/Linux, a daily 50 Hz MAG file is processed in about 25 s.

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
