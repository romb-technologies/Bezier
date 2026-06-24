/*
 * Copyright 2026 Mirko Kokot
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <deque>
#include <iterator>
#include <string>
#include <vector>

#include <nanobind/eigen/dense.h>  // Eigen <-> numpy (Point, Vector, MatrixX2d)
#include <nanobind/nanobind.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/vector.h>  // PointVector, ParamVector, vector<Curve>, ...

#include "Bezier/bezier.h"
#include "Bezier/polycurve.h"
#include "Bezier/utils.h"

namespace nb = nanobind;
using namespace nb::literals;
using namespace Bezier;

// Method names mirror the C++ API (camelCase) so users can cross-reference the C++ docs.
NB_MODULE(bezier, m)
{
  m.doc() = "Python bindings for the Bezier C++ library (romb-technologies). "
            "Points are numpy arrays of shape (2,); point/param lists are (N,2)/(N,).";

  nb::class_<Curve>(m, "Curve")
      .def(nb::init<Eigen::MatrixX2d>(), "points"_a)        // Nx2 control points
      .def(nb::init<const PointVector&>(), "points"_a)      // list of (2,) points
      .def("order", &Curve::order)
      .def("controlPoints", &Curve::controlPoints)
      .def("controlPoint", &Curve::controlPoint, "idx"_a)
      .def("endPoints", &Curve::endPoints)
      .def("polyline", nb::overload_cast<>(&Curve::polyline, nb::const_))
      .def("polyline", nb::overload_cast<double>(&Curve::polyline, nb::const_), "flatness"_a)
      .def("polylineParams", nb::overload_cast<>(&Curve::polylineParams, nb::const_))
      .def("polylineParams", nb::overload_cast<double>(&Curve::polylineParams, nb::const_), "flatness"_a)
      .def("length", nb::overload_cast<>(&Curve::length, nb::const_))
      .def("length", nb::overload_cast<double>(&Curve::length, nb::const_), "t"_a)
      .def("length", nb::overload_cast<double, double>(&Curve::length, nb::const_), "t1"_a, "t2"_a)
      .def("step", &Curve::step, "t"_a, "ds"_a)
      .def("reverse", &Curve::reverse)
      .def("setControlPoint", &Curve::setControlPoint, "idx"_a, "point"_a)
      .def("raiseOrder", &Curve::raiseOrder)
      .def("lowerOrder", &Curve::lowerOrder)
      .def("valueAt", nb::overload_cast<double>(&Curve::valueAt, nb::const_), "t"_a)
      .def("valueAt", nb::overload_cast<const ParamVector&>(&Curve::valueAt, nb::const_), "t_vector"_a)
      .def("curvatureAt", &Curve::curvatureAt, "t"_a)
      .def("curvatureDerivativeAt", &Curve::curvatureDerivativeAt, "t"_a)
      .def("tangentAt", &Curve::tangentAt, "t"_a)
      .def("normalAt", &Curve::normalAt, "t"_a)
      // derivative() returns a const ref into an internal cache; keep the parent alive.
      .def("derivative", nb::overload_cast<>(&Curve::derivative, nb::const_), nb::rv_policy::reference_internal)
      .def("derivative", nb::overload_cast<unsigned>(&Curve::derivative, nb::const_), "n"_a,
           nb::rv_policy::reference_internal)
      .def("derivativeAt", nb::overload_cast<double>(&Curve::derivativeAt, nb::const_), "t"_a)
      .def("derivativeAt", nb::overload_cast<unsigned, double>(&Curve::derivativeAt, nb::const_), "n"_a, "t"_a)
      .def("roots", &Curve::roots)
      .def("extrema", &Curve::extrema)
      // AlignedBox2d has no nanobind caster; return (min, max) instead of a wrapper class.
      .def("boundingBox",
           [](const Curve& c) {
             auto bb = c.boundingBox();
             return std::make_pair(Point(bb.min()), Point(bb.max()));
           })
      .def("splitCurve", nb::overload_cast<const ParamVector&>(&Curve::splitCurve, nb::const_), "t_vector"_a)
      .def("splitCurve", nb::overload_cast<double>(&Curve::splitCurve, nb::const_), nb::arg("t") = 0.5)
      .def("intersections", &Curve::intersections, "curve"_a)
      .def("projectPoint", &Curve::projectPoint, "point"_a)
      .def("distance", &Curve::distance, "point"_a)
      .def("applyContinuity", &Curve::applyContinuity, "curve"_a, "beta_coeffs"_a)
      .def_static("offsetCurve", &Curve::offsetCurve, "curve"_a, "offset"_a, "order"_a = 0)
      .def_static("joinCurves", &Curve::joinCurves, "curve1"_a, "curve2"_a, "order"_a = 0)
      .def_static("fromPolyline", &Curve::fromPolyline, "polyline"_a, "order"_a = 0)
      .def("__repr__", [](const Curve& c) { return "<bezier.Curve order=" + std::to_string(c.order()) + ">"; });

  nb::class_<PolyCurve>(m, "PolyCurve")
      .def(nb::init<>())
      // PolyCurve takes std::deque<Curve> (no deque caster); accept a list and build the deque.
      .def(
          "__init__",
          [](PolyCurve* self, std::vector<Curve> curves) {
            new (self) PolyCurve(
                std::deque<Curve>(std::make_move_iterator(curves.begin()), std::make_move_iterator(curves.end())));
          },
          "curves"_a)
      .def("insertAt", &PolyCurve::insertAt, "idx"_a, "curve"_a)
      .def("insertFront", &PolyCurve::insertFront, "curve"_a)
      .def("insertBack", &PolyCurve::insertBack, "curve"_a)
      .def("removeAt", &PolyCurve::removeAt, "idx"_a)
      .def("removeFirst", &PolyCurve::removeFirst)
      .def("removeBack", &PolyCurve::removeBack)
      .def("size", &PolyCurve::size)
      .def("__len__", &PolyCurve::size)
      .def("curveIdx", &PolyCurve::curveIdx, "t"_a)
      // curve(idx) returns a reference into the internal deque; keep the parent alive.
      .def("curve", nb::overload_cast<unsigned>(&PolyCurve::curve), "idx"_a, nb::rv_policy::reference_internal)
      .def("polyline", nb::overload_cast<>(&PolyCurve::polyline, nb::const_))
      .def("polyline", nb::overload_cast<double>(&PolyCurve::polyline, nb::const_), "flatness"_a)
      .def("polylineParams", nb::overload_cast<>(&PolyCurve::polylineParams, nb::const_))
      .def("polylineParams", nb::overload_cast<double>(&PolyCurve::polylineParams, nb::const_), "flatness"_a)
      .def("length", nb::overload_cast<>(&PolyCurve::length, nb::const_))
      .def("length", nb::overload_cast<double>(&PolyCurve::length, nb::const_), "t"_a)
      .def("length", nb::overload_cast<double, double>(&PolyCurve::length, nb::const_), "t1"_a, "t2"_a)
      .def("step", &PolyCurve::step, "t"_a, "ds"_a)
      .def("endPoints", &PolyCurve::endPoints)
      .def("controlPoints", &PolyCurve::controlPoints)
      .def("setControlPoint", &PolyCurve::setControlPoint, "idx"_a, "point"_a)
      .def("valueAt", nb::overload_cast<double>(&PolyCurve::valueAt, nb::const_), "t"_a)
      .def("valueAt", nb::overload_cast<const ParamVector&>(&PolyCurve::valueAt, nb::const_), "t_vector"_a)
      .def("curvatureAt", &PolyCurve::curvatureAt, "t"_a)
      .def("curvatureDerivativeAt", &PolyCurve::curvatureDerivativeAt, "t"_a)
      .def("tangentAt", &PolyCurve::tangentAt, "t"_a)
      .def("normalAt", &PolyCurve::normalAt, "t"_a)
      .def("derivativeAt", nb::overload_cast<double>(&PolyCurve::derivativeAt, nb::const_), "t"_a)
      .def("derivativeAt", nb::overload_cast<unsigned, double>(&PolyCurve::derivativeAt, nb::const_), "n"_a, "t"_a)
      .def("boundingBox",
           [](const PolyCurve& c) {
             auto bb = c.boundingBox();
             return std::make_pair(Point(bb.min()), Point(bb.max()));
           })
      .def("intersections", nb::overload_cast<const Curve&>(&PolyCurve::intersections, nb::const_), "curve"_a)
      .def("intersections", nb::overload_cast<const PolyCurve&>(&PolyCurve::intersections, nb::const_), "poly_curve"_a)
      .def("projectPoint", nb::overload_cast<const Point&>(&PolyCurve::projectPoint, nb::const_), "point"_a)
      .def("projectPoint", nb::overload_cast<const PointVector&>(&PolyCurve::projectPoint, nb::const_),
           "point_vector"_a)
      .def("distance", nb::overload_cast<const Point&>(&PolyCurve::distance, nb::const_), "point"_a)
      .def("distance", nb::overload_cast<const PointVector&>(&PolyCurve::distance, nb::const_), "point_vector"_a)
      .def("__repr__",
           [](const PolyCurve& c) { return "<bezier.PolyCurve size=" + std::to_string(c.size()) + ">"; });

  // only the non-inline Utils symbols; inline math helpers stay internal.
  auto utils = m.def_submodule("Utils");
  utils.def("visvalingamWyatt", &Utils::visvalingamWyatt, "polyline"_a);
  utils.def("solvePolynomial", &Utils::solvePolynomial, "polynomial"_a);
  utils.def("fitBezier", &Utils::fitBezier, "points"_a, "order"_a);
}
