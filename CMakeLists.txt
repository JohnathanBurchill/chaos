# CHAOS: CMakeLists.txt

# Copyright (C) 2023  Johnathan K Burchill

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

CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

PROJECT(chaos)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMake")

if(CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Debugging enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -std=gnu99 ")
else (CMAKE_BUILD_TYPE STREQUAL Debug)
    message( "-- Optimization enabled")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -std=gnu99 ")
endif(CMAKE_BUILD_TYPE STREQUAL Debug)


# requires -lm on linux
FIND_LIBRARY(MATH m)
FIND_LIBRARY(CDF cdf)
FIND_LIBRARY(GSL gsl)
FIND_LIBRARY(CBLAS gslcblas)


INCLUDE_DIRECTORIES(include)

ADD_LIBRARY(chaostrace trace.c model.c shc.c)

ADD_EXECUTABLE(chaos chaos.c cdf_utils.c cdf_vars.c cdf_attrs.c shc.c model.c)
TARGET_LINK_LIBRARIES(chaos ${LIBS} ${CDF} -lgsl -lm -lgslcblas)

ADD_EXECUTABLE(tracechaos tracechaos.c)
TARGET_LINK_LIBRARIES(tracechaos chaostrace ${LIBS} -lgsl -lm -lgslcblas)

ADD_EXECUTABLE(themis_asi_fieldlines themis_asi_fieldlines.c util.c)
TARGET_LINK_LIBRARIES(themis_asi_fieldlines chaostrace ${LIBS} ${CDF} -lgsl -lm -lgslcblas)

ADD_EXECUTABLE(chaos_calc chaos_calc.c util.c)
TARGET_LINK_LIBRARIES(chaos_calc chaostrace ${LIBS} -lgsl -lm -lgslcblas)

install(TARGETS chaos DESTINATION $ENV{HOME}/bin)
install(TARGETS tracechaos DESTINATION $ENV{HOME}/bin)
install(TARGETS themis_asi_fieldlines DESTINATION $ENV{HOME}/bin)
install(TARGETS chaos_calc DESTINATION $ENV{HOME}/bin)
