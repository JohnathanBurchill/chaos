# CHAOS: Docker/Dockerfile

# Copyright (C) 2022  Johnathan K Burchill

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# adapted from TII Ion Drift Processor

FROM ubuntu
RUN apt -y update
RUN apt -y upgrade
RUN apt -y install wget
RUN apt -y install gcc
RUN apt -y install g++
RUN apt -y install cmake
RUN apt -y install make
RUN apt -y install bash
RUN apt -y install zlib1g
RUN apt -y install zlib1g-dev
RUN apt -y install git
RUN apt -y install gdb
RUN apt -y install cgdb
RUN apt -y install libgsl-dev
RUN apt -y install vim
RUN apt -y install pkg-config
RUN addgroup science
RUN useradd -G science dataflow
USER dataflow
WORKDIR /home/dataflow/src
RUN cd /home/dataflow && \
wget https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf38_1/linux/cdf38_1-dist-cdf.tar.gz && \
tar xf cdf38_1-dist-cdf.tar.gz && \
cd cdf38_1-dist && \
make OS=linux ENV=gnu CURSES=no all
USER root
RUN cd /home/dataflow/cdf38_1-dist && \
make INSTALLDIR=/usr/local install && \
ldconfig
USER dataflow
RUN cd /home/dataflow && \
rm -rf cdf38_1-dist*
VOLUME ["/home/dataflow/src"]
ENTRYPOINT ["/bin/bash"]



