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

#include "bezier_c.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <string>
#include <vector>

#include "Bezier/bezier.h"
#include "Bezier/polycurve.h"
#include "Bezier/utils.h"

using Bezier::Curve;
using Bezier::Point;
using Bezier::PointVector;
using Bezier::PolyCurve;

namespace
{
thread_local std::string g_last_error;

// Run f() inside the exception boundary; on throw, record the message and return a default.
template <class F> auto guard(F&& f) -> decltype(f())
{
  try
  {
    g_last_error.clear();
    return f();
  }
  catch (const std::exception& e)
  {
    g_last_error = e.what();
  }
  catch (...)
  {
    g_last_error = "unknown error";
  }
  return decltype(f())();
}

Curve* cv(BzCurve* h) { return reinterpret_cast<Curve*>(h); }
const Curve* cv(const BzCurve* h) { return reinterpret_cast<const Curve*>(h); }
BzCurve* hc(Curve* c) { return reinterpret_cast<BzCurve*>(c); }

PolyCurve* pv(BzPolyCurve* h) { return reinterpret_cast<PolyCurve*>(h); }
const PolyCurve* pv(const BzPolyCurve* h) { return reinterpret_cast<const PolyCurve*>(h); }
BzPolyCurve* hp(PolyCurve* p) { return reinterpret_cast<BzPolyCurve*>(p); }

PointVector toPoints(const double* xy, int n)
{
  PointVector pts;
  pts.reserve(n);
  for (int k = 0; k < n; ++k)
    pts.emplace_back(xy[2 * k], xy[2 * k + 1]);
  return pts;
}

template <class T> T* allocCopy(const std::vector<T>& v, int* out_count)
{
  *out_count = static_cast<int>(v.size());
  if (v.empty())
    return nullptr;
  T* buf = static_cast<T*>(std::malloc(sizeof(T) * v.size()));
  std::copy(v.begin(), v.end(), buf);
  return buf;
}

double* allocPoints(const PointVector& pts, int* out_count)
{
  *out_count = static_cast<int>(pts.size());
  if (pts.empty())
    return nullptr;
  double* buf = static_cast<double*>(std::malloc(sizeof(double) * 2 * pts.size()));
  for (size_t k = 0; k < pts.size(); ++k)
  {
    buf[2 * k] = pts[k].x();
    buf[2 * k + 1] = pts[k].y();
  }
  return buf;
}

void writePoint(const Point& p, double out[2])
{
  out[0] = p.x();
  out[1] = p.y();
}
} // namespace

extern "C"
{

const char* bz_last_error(void) { return g_last_error.empty() ? nullptr : g_last_error.c_str(); }
void bz_clear_error(void) { g_last_error.clear(); }
void bz_free(void* p) { std::free(p); }

/* ============================ Curve ============================ */

BzCurve* bz_curve_new(const double* xy, int n_points)
{
  return guard([&] { return hc(new Curve(toPoints(xy, n_points))); });
}
BzCurve* bz_curve_copy(const BzCurve* c) { return guard([&] { return hc(new Curve(*cv(c))); }); }
void bz_curve_destroy(BzCurve* c) { delete cv(c); }

unsigned bz_curve_order(const BzCurve* c) { return guard([&] { return cv(c)->order(); }); }

double* bz_curve_control_points(const BzCurve* c, int* out_count)
{
  return guard([&] { return allocPoints(cv(c)->controlPoints(), out_count); });
}
void bz_curve_control_point(const BzCurve* c, unsigned idx, double out[2])
{
  guard([&] { writePoint(cv(c)->controlPoint(idx), out); });
}
void bz_curve_end_points(const BzCurve* c, double out_first[2], double out_second[2])
{
  guard([&] {
    auto ep = cv(c)->endPoints();
    writePoint(ep.first, out_first);
    writePoint(ep.second, out_second);
  });
}

double* bz_curve_polyline(const BzCurve* c, int use_flatness, double flatness, int* out_count)
{
  return guard([&] { return allocPoints(use_flatness ? cv(c)->polyline(flatness) : cv(c)->polyline(), out_count); });
}
double* bz_curve_polyline_params(const BzCurve* c, int use_flatness, double flatness, int* out_count)
{
  return guard(
      [&] { return allocCopy(use_flatness ? cv(c)->polylineParams(flatness) : cv(c)->polylineParams(), out_count); });
}

double bz_curve_length(const BzCurve* c) { return guard([&] { return cv(c)->length(); }); }
double bz_curve_length_to(const BzCurve* c, double t) { return guard([&] { return cv(c)->length(t); }); }
double bz_curve_length_between(const BzCurve* c, double t1, double t2)
{
  return guard([&] { return cv(c)->length(t1, t2); });
}
double bz_curve_step(const BzCurve* c, double t, double ds) { return guard([&] { return cv(c)->step(t, ds); }); }

void bz_curve_reverse(BzCurve* c) { guard([&] { cv(c)->reverse(); }); }
void bz_curve_set_control_point(BzCurve* c, unsigned idx, const double point[2])
{
  guard([&] { cv(c)->setControlPoint(idx, Point(point[0], point[1])); });
}
void bz_curve_raise_order(BzCurve* c) { guard([&] { cv(c)->raiseOrder(); }); }
void bz_curve_lower_order(BzCurve* c) { guard([&] { cv(c)->lowerOrder(); }); }

void bz_curve_value_at(const BzCurve* c, double t, double out[2])
{
  guard([&] { writePoint(cv(c)->valueAt(t), out); });
}
double* bz_curve_value_at_many(const BzCurve* c, const double* t, int n, int* out_count)
{
  return guard([&] {
    Bezier::ParamVector tv(t, t + n);
    Eigen::MatrixX2d m = cv(c)->valueAt(tv);
    PointVector pts(m.rows());
    for (Eigen::Index r = 0; r < m.rows(); ++r)
      pts[r] = m.row(r).transpose();
    return allocPoints(pts, out_count);
  });
}
double bz_curve_curvature_at(const BzCurve* c, double t) { return guard([&] { return cv(c)->curvatureAt(t); }); }
double bz_curve_curvature_derivative_at(const BzCurve* c, double t)
{
  return guard([&] { return cv(c)->curvatureDerivativeAt(t); });
}
void bz_curve_tangent_at(const BzCurve* c, double t, double out[2])
{
  guard([&] { writePoint(cv(c)->tangentAt(t), out); });
}
void bz_curve_normal_at(const BzCurve* c, double t, double out[2])
{
  guard([&] { writePoint(cv(c)->normalAt(t), out); });
}

BzCurve* bz_curve_derivative(const BzCurve* c, unsigned n)
{
  // copy out of the parent's cache so the C# handle owns an independent curve
  return guard([&] { return hc(new Curve(cv(c)->derivative(n))); });
}
void bz_curve_derivative_at(const BzCurve* c, unsigned n, double t, double out[2])
{
  guard([&] { writePoint(cv(c)->derivativeAt(n, t), out); });
}

double* bz_curve_roots(const BzCurve* c, int* out_count)
{
  return guard([&] { return allocCopy(cv(c)->roots(), out_count); });
}
double* bz_curve_extrema(const BzCurve* c, int* out_count)
{
  return guard([&] { return allocCopy(cv(c)->extrema(), out_count); });
}
void bz_curve_bounding_box(const BzCurve* c, double out_min[2], double out_max[2])
{
  guard([&] {
    auto bb = cv(c)->boundingBox();
    writePoint(bb.min(), out_min);
    writePoint(bb.max(), out_max);
  });
}

BzCurve** bz_curve_split(const BzCurve* c, const double* t, int n, int* out_count)
{
  return guard([&]() -> BzCurve** {
    Bezier::ParamVector tv(t, t + n);
    std::vector<Curve> parts = cv(c)->splitCurve(tv);
    *out_count = static_cast<int>(parts.size());
    if (parts.empty())
      return nullptr;
    BzCurve** buf = static_cast<BzCurve**>(std::malloc(sizeof(BzCurve*) * parts.size()));
    for (size_t k = 0; k < parts.size(); ++k)
      buf[k] = hc(new Curve(std::move(parts[k])));
    return buf;
  });
}
void bz_curve_split_at(const BzCurve* c, double t, BzCurve** out_left, BzCurve** out_right)
{
  guard([&] {
    auto pr = cv(c)->splitCurve(t);
    *out_left = hc(new Curve(std::move(pr.first)));
    *out_right = hc(new Curve(std::move(pr.second)));
  });
}

double* bz_curve_intersections(const BzCurve* c, const BzCurve* other, int* out_count)
{
  return guard([&] { return allocPoints(cv(c)->intersections(*cv(other)), out_count); });
}
double bz_curve_project_point(const BzCurve* c, const double point[2])
{
  return guard([&] { return cv(c)->projectPoint(Point(point[0], point[1])); });
}
double bz_curve_distance(const BzCurve* c, const double point[2])
{
  return guard([&] { return cv(c)->distance(Point(point[0], point[1])); });
}
void bz_curve_apply_continuity(BzCurve* c, const BzCurve* other, const double* beta, int n_beta)
{
  guard([&] { cv(c)->applyContinuity(*cv(other), std::vector<double>(beta, beta + n_beta)); });
}

BzCurve* bz_curve_offset(const BzCurve* c, double offset, unsigned order)
{
  return guard([&] { return hc(new Curve(Curve::offsetCurve(*cv(c), offset, order))); });
}
BzCurve* bz_curve_join(const BzCurve* c1, const BzCurve* c2, unsigned order)
{
  return guard([&] { return hc(new Curve(Curve::joinCurves(*cv(c1), *cv(c2), order))); });
}
BzCurve* bz_curve_from_polyline(const double* xy, int n_points, unsigned order)
{
  return guard([&] { return hc(new Curve(Curve::fromPolyline(toPoints(xy, n_points), order))); });
}

/* ============================ PolyCurve ============================ */

BzPolyCurve* bz_polycurve_new(void) { return guard([&] { return hp(new PolyCurve()); }); }
BzPolyCurve* bz_polycurve_new_from(BzCurve* const* curves, int n)
{
  return guard([&] {
    std::deque<Curve> dq;
    for (int k = 0; k < n; ++k)
      dq.push_back(*cv(curves[k]));
    return hp(new PolyCurve(std::move(dq)));
  });
}
BzPolyCurve* bz_polycurve_copy(const BzPolyCurve* p) { return guard([&] { return hp(new PolyCurve(*pv(p))); }); }
void bz_polycurve_destroy(BzPolyCurve* p) { delete pv(p); }

void bz_polycurve_insert_at(BzPolyCurve* p, unsigned idx, const BzCurve* c)
{
  guard([&] { pv(p)->insertAt(idx, *cv(c)); });
}
void bz_polycurve_insert_front(BzPolyCurve* p, const BzCurve* c) { guard([&] { pv(p)->insertFront(*cv(c)); }); }
void bz_polycurve_insert_back(BzPolyCurve* p, const BzCurve* c) { guard([&] { pv(p)->insertBack(*cv(c)); }); }
void bz_polycurve_remove_at(BzPolyCurve* p, unsigned idx) { guard([&] { pv(p)->removeAt(idx); }); }
void bz_polycurve_remove_first(BzPolyCurve* p) { guard([&] { pv(p)->removeFirst(); }); }
void bz_polycurve_remove_back(BzPolyCurve* p) { guard([&] { pv(p)->removeBack(); }); }

unsigned bz_polycurve_size(const BzPolyCurve* p) { return guard([&] { return pv(p)->size(); }); }
unsigned bz_polycurve_curve_idx(const BzPolyCurve* p, double t) { return guard([&] { return pv(p)->curveIdx(t); }); }
BzCurve* bz_polycurve_curve(const BzPolyCurve* p, unsigned idx)
{
  return guard([&] { return hc(new Curve(pv(p)->curve(idx))); });
}

double* bz_polycurve_polyline(const BzPolyCurve* p, int use_flatness, double flatness, int* out_count)
{
  return guard([&] { return allocPoints(use_flatness ? pv(p)->polyline(flatness) : pv(p)->polyline(), out_count); });
}
double* bz_polycurve_polyline_params(const BzPolyCurve* p, int use_flatness, double flatness, int* out_count)
{
  return guard(
      [&] { return allocCopy(use_flatness ? pv(p)->polylineParams(flatness) : pv(p)->polylineParams(), out_count); });
}

double bz_polycurve_length(const BzPolyCurve* p) { return guard([&] { return pv(p)->length(); }); }
double bz_polycurve_length_to(const BzPolyCurve* p, double t) { return guard([&] { return pv(p)->length(t); }); }
double bz_polycurve_length_between(const BzPolyCurve* p, double t1, double t2)
{
  return guard([&] { return pv(p)->length(t1, t2); });
}
double bz_polycurve_step(const BzPolyCurve* p, double t, double ds) { return guard([&] { return pv(p)->step(t, ds); }); }

void bz_polycurve_end_points(const BzPolyCurve* p, double out_first[2], double out_second[2])
{
  guard([&] {
    auto ep = pv(p)->endPoints();
    writePoint(ep.first, out_first);
    writePoint(ep.second, out_second);
  });
}
double* bz_polycurve_control_points(const BzPolyCurve* p, int* out_count)
{
  return guard([&] { return allocPoints(pv(p)->controlPoints(), out_count); });
}
void bz_polycurve_set_control_point(BzPolyCurve* p, unsigned idx, const double point[2])
{
  guard([&] { pv(p)->setControlPoint(idx, Point(point[0], point[1])); });
}

void bz_polycurve_value_at(const BzPolyCurve* p, double t, double out[2])
{
  guard([&] { writePoint(pv(p)->valueAt(t), out); });
}
double* bz_polycurve_value_at_many(const BzPolyCurve* p, const double* t, int n, int* out_count)
{
  return guard([&] {
    Bezier::ParamVector tv(t, t + n);
    return allocPoints(pv(p)->valueAt(tv), out_count);
  });
}
double bz_polycurve_curvature_at(const BzPolyCurve* p, double t) { return guard([&] { return pv(p)->curvatureAt(t); }); }
double bz_polycurve_curvature_derivative_at(const BzPolyCurve* p, double t)
{
  return guard([&] { return pv(p)->curvatureDerivativeAt(t); });
}
void bz_polycurve_tangent_at(const BzPolyCurve* p, double t, double out[2])
{
  guard([&] { writePoint(pv(p)->tangentAt(t), out); });
}
void bz_polycurve_normal_at(const BzPolyCurve* p, double t, double out[2])
{
  guard([&] { writePoint(pv(p)->normalAt(t), out); });
}
void bz_polycurve_derivative_at(const BzPolyCurve* p, unsigned n, double t, double out[2])
{
  guard([&] { writePoint(pv(p)->derivativeAt(n, t), out); });
}
void bz_polycurve_bounding_box(const BzPolyCurve* p, double out_min[2], double out_max[2])
{
  guard([&] {
    auto bb = pv(p)->boundingBox();
    writePoint(bb.min(), out_min);
    writePoint(bb.max(), out_max);
  });
}

double* bz_polycurve_intersections_curve(const BzPolyCurve* p, const BzCurve* c, int* out_count)
{
  return guard([&] { return allocPoints(pv(p)->intersections(*cv(c)), out_count); });
}
double* bz_polycurve_intersections_poly(const BzPolyCurve* p, const BzPolyCurve* other, int* out_count)
{
  return guard([&] { return allocPoints(pv(p)->intersections(*pv(other)), out_count); });
}
double bz_polycurve_project_point(const BzPolyCurve* p, const double point[2])
{
  return guard([&] { return pv(p)->projectPoint(Point(point[0], point[1])); });
}
double* bz_polycurve_project_points(const BzPolyCurve* p, const double* xy, int n_points, int* out_count)
{
  return guard([&] { return allocCopy(pv(p)->projectPoint(toPoints(xy, n_points)), out_count); });
}
double bz_polycurve_distance(const BzPolyCurve* p, const double point[2])
{
  return guard([&] { return pv(p)->distance(Point(point[0], point[1])); });
}
double* bz_polycurve_distances(const BzPolyCurve* p, const double* xy, int n_points, int* out_count)
{
  return guard([&] { return allocCopy(pv(p)->distance(toPoints(xy, n_points)), out_count); });
}

/* ============================ Utils ============================ */

int* bz_utils_visvalingam_wyatt(const double* xy, int n_points, int* out_count)
{
  return guard([&]() -> int* {
    std::vector<unsigned> idx = Bezier::Utils::visvalingamWyatt(toPoints(xy, n_points));
    *out_count = static_cast<int>(idx.size());
    if (idx.empty())
      return nullptr;
    int* buf = static_cast<int*>(std::malloc(sizeof(int) * idx.size()));
    for (size_t k = 0; k < idx.size(); ++k)
      buf[k] = static_cast<int>(idx[k]);
    return buf;
  });
}
double* bz_utils_solve_polynomial(const double* coeffs, int n, int* out_count)
{
  return guard([&] {
    Eigen::VectorXd poly = Eigen::Map<const Eigen::VectorXd>(coeffs, n);
    return allocCopy(Bezier::Utils::solvePolynomial(poly), out_count);
  });
}
BzCurve* bz_utils_fit_bezier(const double* xy, int n_points, unsigned order)
{
  return guard([&] { return hc(new Curve(Bezier::Utils::fitBezier(toPoints(xy, n_points), order))); });
}

} // extern "C"
