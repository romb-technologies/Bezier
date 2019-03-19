# Bezier-Cpp
[![Build Status](https://travis-ci.com/stribor14/Bezier-cpp.svg?branch=master)](https://travis-ci.com/stribor14/Bezier-cpp)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/aceb46ce7de1407abd56cfc127dba5f1)](https://www.codacy.com/app/stribor14/Bezier-cpp?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=stribor14/Bezier-cpp&amp;utm_campaign=Badge_Grade)

Fast and lightweight class for using the Bezier curves of any order in C++

*Algorithm implementations are based on [A Primer on Bezier Curves](https://pomax.github.io/bezierinfo/) by Pomax*

## Key Features
  - dynamic operations
  - any order of curve
  - optimized for real-time calculations

## Implementeed methods
  - get value, curvature, tangent and normal for parameter *t*
  - get t from projection any point onto a curve
  - get derivative curve
  - split into two subcurves
  - find extremes and bounding box
  - find points of intersection with another curve
  - elevate/lower order
  - manipulate control points
  - manipulate dot on curve (only for quadratic and cubic curves)

## Dependencies
  - c++11
  - Eigen3

*Add compile flag* `-march=native` *when compiling to use vectorization with Eigen.*

## Example program
A small Qt5 based program written as a playground for manipulating Bezier curves.
### Usage
 - starts with two Bezier curves (with 4 and 5 control points respectively)
 - Zoom in/out: *__Ctrl + mouse wheel__*
 - Manipulate control point or point on curve: *__Left mouse buttom__*
 - Project mouse pointer on all curves and show tangent: *__Right mouse buttom__*
 - Split curve at mouse point: *__Middle mouse buttom__*
 - Raise order of the curve: *__Double left click__*
 - Lower order of the curve *__Double right click__*
 - Toggle bounding boxes and curve intersections: *__Double middle click__*

### Additional dependencies
 - qt5-default 

## Licence
Apache License Version 2.0
