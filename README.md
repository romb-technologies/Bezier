# Bezier-Cpp
[![Build Status](https://travis-ci.com/stribor14/Bezier-cpp.svg?branch=master)](https://travis-ci.com/stribor14/Bezier-cpp)
<<<<<<< HEAD
![v0.1](https://img.shields.io/badge/version-0.1-blue.svg)
=======
![v0.2](https://img.shields.io/badge/version-0.2-blue.svg)
>>>>>>> pre_release
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/aceb46ce7de1407abd56cfc127dba5f1)](https://www.codacy.com/app/stribor14/Bezier-cpp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=stribor14/Bezier-cpp&amp;utm_campaign=Badge_Grade)

Fast and lightweight class for using the Bezier curves of any order in C++

*Algorithm implementations are based on [A Primer on Bezier Curves](https://pomax.github.io/bezierinfo/) by Pomax*

## Key Features
  - Any number of control points
  - Fast operations on curves
<<<<<<< HEAD
  - Dynamic manipulation  

## Implemented methods
  - Get value, curvature, tangent and normal for parameter *t*
  - Get t from projection any point onto a curve
  - Get derivative curve
  - Split into two subcurves
  - Find extremes and bounding box
  - Find points of intersection with another curve
  - Elevate/lower order
  - Manipulate control points
  - Manipulate dot on curve (only for quadratic and cubic curves)
  
## In development
  - <img src="https://img.shields.io/badge/v.0.2-indev-yellow.svg" alt="v0.2 indev" align="top"> Bezier polycurves
    - [x] Create polycurve
    - [x] Equivalent methods from Bezier curves
    - [x] Adding and removing curves
    - [x] Curve continuity
      - [x] Parametric continuity
      - [x] Geometric continuity
    - [x] Document polycurve features
    - [ ] More sophisticated example
  - <img src="https://img.shields.io/badge/v.0.3-planned-red.svg" alt="v0.3 planned" align="top"> Bezier shapes
=======
  - Dynamic manipulation
  - Composite Bezier curves (polycurves)

CMake *find_package()* compatible!
```
find_package(Bezier)
target_link_libraries(target bezier)
```

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
  
## Wish list
    - [ ] Polycurve - oversee continuities between consecutive sub-curves
    - [ ] Polycurve - propagation of sub-curve manipulation depending on continutiy
    - [ ] Bezier shapes
    - [ ] More sophisticated example
>>>>>>> pre_release

## Dependencies
  - c++11
  - Eigen3

## Instalation
### System-wide installation
```
git clone https://github.com/stribor14/Bezier-cpp
mkdir Bezier-cpp/build
cd Bezier-cpp/build
cmake ..
make
make install
```
### ROS
- for use within a ROS workspace without the system-wide installation, clone the repo to src folder in you catkin workspace 

## Example program __[OUTDATED]__ 
A small Qt5 based program written as a playground for manipulating Bezier curves.
<<<<<<< HEAD
### Usage
 - Starts with two Bezier curves (with 4 and 5 control points respectively)
 - Zoom in/out: *__Ctrl + mouse wheel__*
 - Manipulate control point or point on curve: *__Left mouse buttom__*
 - Project mouse pointer on all curves and show tangent: *__Right mouse buttom__*
 - Split curve at mouse point: *__Middle mouse buttom__*
 - Raise order of the curve: *__Double left click__*
 - Lower order of the curve *__Double right click__*
 - Toggle bounding boxes and curve intersections: *__Double middle click__*
=======
- press *__H__* for a list of possible actions
>>>>>>> pre_release

### Additional dependencies
 - qt5-default 

## Licence
Apache License Version 2.0
