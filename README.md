# Bezier-Cpp
[![Build Status](https://travis-ci.com/stribor14/Bezier-cpp.svg?branch=master)](https://travis-ci.com/stribor14/Bezier-cpp)
![v0.2.0](https://img.shields.io/badge/version-0.2.0-blue.svg)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/aceb46ce7de1407abd56cfc127dba5f1)](https://www.codacy.com/app/stribor14/Bezier-cpp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=stribor14/Bezier-cpp&amp;utm_campaign=Badge_Grade)

Fast and lightweight class for using the Bezier curves of any order in C++

*Algorithm implementations are based on [A Primer on Bezier Curves](https://pomax.github.io/bezierinfo/) by Pomax*

## Key Features
  - Any number of control points
  - Fast operations on curves
  - Dynamic manipulation
  - Composite Bezier curves (polycurves)

## Implemented methods
  - Get value, derivative, curvature, tangent and normal for parameter *t*
  - Get t from projection any point onto a curve
  - Get precise length for any part of curve
  - Get a derivative curve (hodograph)
  - Split into two subcurves
  - Find curve roots and bounding box
  - Find points of intersection with another curve
  - Elevate/lower order
  - Apply parametric and geometric continuities
  - etc.
  
## In development
  - <img src="https://img.shields.io/badge/v.0.2.1-planned-yellow.svg" alt="v0.2.1 indev" align="top"> Bezier polycurves
    - [ ] Polycurve - oversee continuities between consecutive sub-curves
    - [ ] Polycurve - propagation of sub-curve manipulation depending on continutiy
    - [ ] More sophisticated example
  - <img src="https://img.shields.io/badge/v.0.3-planned-red.svg" alt="v0.3 planned" align="top"> Bezier shapes

## Dependencies
  - c++11
  - Eigen3

## Instalation
CMake *find_package()* compatible!
```
find_package(BezierCpp)
target_link_libraries(target BezierCpp)
```
### System-wide installation
```
git clone https://github.com/stribor14/Bezier-cpp
mkdir Bezicer-cpp/build
cd Bezicer-cpp/build
cmake ..
make
make install
```
### ROS
- for use within a ROS workspace without the system-wide installation, clone the repo to src folder in you catkin workspace 

## Example program __[OUTDATED]__ 
A small Qt5 based program written as a playground for manipulating Bezier curves.
- press *__H__* for a list of possible actions

### Additional dependencies
 - qt5-default 

## Licence
Apache License Version 2.0
