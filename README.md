# Bezier-Cpp
Fast and lightweight class for using the Bezier curves of any order in C++
# Key Features:
  - dynamic operations
  - any order of curve
  - optimized for real-time calculations

# Implementeed methods:
  - get value, curvature, tangent and normal for parameter *t*
  - get t from projection any point onto a curve
  - get derivative curve
  - split into two subcurves
  - find extremes and bounding box
  - find points of intersection with another curve
  - elevate/lower order
  - manipulate control points
  - manipulate dot on curve

# Dependencies:
- c++11
- Eigen3

*Add compile flag* `-march=native` *when compiling to use vectorization with Eigen.*

# Licence
Apache License Version 2.0
