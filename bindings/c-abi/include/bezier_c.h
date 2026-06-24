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

/*
 * Flat C ABI over Bezier::Curve / Bezier::PolyCurve / Bezier::Utils, for P/Invoke.
 *
 * Conventions:
 *  - Handles (BzCurve*, BzPolyCurve*) are owned by the caller; release with the matching
 *    *_destroy. Functions that return a handle transfer ownership to the caller.
 *  - A single point is written into a caller-supplied double[2] (x, y).
 *  - A variable-length list is returned as a malloc'd buffer + *out_count, freed by the
 *    caller with bz_free. Point lists are interleaved x,y (2*count doubles); param/scalar
 *    lists are count doubles. A curve list is a malloc'd array of count BzCurve* (free the
 *    array with bz_free; each curve is owned by the caller).
 *  - C++ exceptions never cross this boundary: on failure a function returns a sentinel
 *    (null / 0) and sets bz_last_error(). Callers should clear/check it around fallible calls.
 */

#ifndef BEZIER_C_H
#define BEZIER_C_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct BzCurve BzCurve;
typedef struct BzPolyCurve BzPolyCurve;

/* --- error handling --- */
const char* bz_last_error(void);          /* null if no error since last clear */
void bz_clear_error(void);
void bz_free(void* p);                    /* frees buffers returned by the functions below */

/* ============================ Curve ============================ */

BzCurve* bz_curve_new(const double* xy, int n_points); /* xy: interleaved, 2*n_points doubles */
BzCurve* bz_curve_copy(const BzCurve* c);
void bz_curve_destroy(BzCurve* c);

unsigned bz_curve_order(const BzCurve* c);
double* bz_curve_control_points(const BzCurve* c, int* out_count);     /* points */
void bz_curve_control_point(const BzCurve* c, unsigned idx, double out[2]);
void bz_curve_end_points(const BzCurve* c, double out_first[2], double out_second[2]);

double* bz_curve_polyline(const BzCurve* c, int use_flatness, double flatness, int* out_count);
double* bz_curve_polyline_params(const BzCurve* c, int use_flatness, double flatness, int* out_count);

double bz_curve_length(const BzCurve* c);
double bz_curve_length_to(const BzCurve* c, double t);
double bz_curve_length_between(const BzCurve* c, double t1, double t2);
double bz_curve_step(const BzCurve* c, double t, double ds);

void bz_curve_reverse(BzCurve* c);
void bz_curve_set_control_point(BzCurve* c, unsigned idx, const double point[2]);
void bz_curve_raise_order(BzCurve* c);
void bz_curve_lower_order(BzCurve* c);

void bz_curve_value_at(const BzCurve* c, double t, double out[2]);
double* bz_curve_value_at_many(const BzCurve* c, const double* t, int n, int* out_count); /* points */
double bz_curve_curvature_at(const BzCurve* c, double t);
double bz_curve_curvature_derivative_at(const BzCurve* c, double t);
void bz_curve_tangent_at(const BzCurve* c, double t, double out[2]);
void bz_curve_normal_at(const BzCurve* c, double t, double out[2]);

BzCurve* bz_curve_derivative(const BzCurve* c, unsigned n); /* n>=1; returns an owned copy */
void bz_curve_derivative_at(const BzCurve* c, unsigned n, double t, double out[2]); /* n>=1 */

double* bz_curve_roots(const BzCurve* c, int* out_count);
double* bz_curve_extrema(const BzCurve* c, int* out_count);
void bz_curve_bounding_box(const BzCurve* c, double out_min[2], double out_max[2]);

BzCurve** bz_curve_split(const BzCurve* c, const double* t, int n, int* out_count);
void bz_curve_split_at(const BzCurve* c, double t, BzCurve** out_left, BzCurve** out_right);

double* bz_curve_intersections(const BzCurve* c, const BzCurve* other, int* out_count); /* points */
double bz_curve_project_point(const BzCurve* c, const double point[2]);
double bz_curve_distance(const BzCurve* c, const double point[2]);
void bz_curve_apply_continuity(BzCurve* c, const BzCurve* other, const double* beta, int n_beta);

BzCurve* bz_curve_offset(const BzCurve* c, double offset, unsigned order);
BzCurve* bz_curve_join(const BzCurve* c1, const BzCurve* c2, unsigned order);
BzCurve* bz_curve_from_polyline(const double* xy, int n_points, unsigned order);

/* ============================ PolyCurve ============================ */

BzPolyCurve* bz_polycurve_new(void);
BzPolyCurve* bz_polycurve_new_from(BzCurve* const* curves, int n); /* copies each curve */
BzPolyCurve* bz_polycurve_copy(const BzPolyCurve* p);
void bz_polycurve_destroy(BzPolyCurve* p);

void bz_polycurve_insert_at(BzPolyCurve* p, unsigned idx, const BzCurve* c);
void bz_polycurve_insert_front(BzPolyCurve* p, const BzCurve* c);
void bz_polycurve_insert_back(BzPolyCurve* p, const BzCurve* c);
void bz_polycurve_remove_at(BzPolyCurve* p, unsigned idx);
void bz_polycurve_remove_first(BzPolyCurve* p);
void bz_polycurve_remove_back(BzPolyCurve* p);

unsigned bz_polycurve_size(const BzPolyCurve* p);
unsigned bz_polycurve_curve_idx(const BzPolyCurve* p, double t);
BzCurve* bz_polycurve_curve(const BzPolyCurve* p, unsigned idx); /* owned copy */

double* bz_polycurve_polyline(const BzPolyCurve* p, int use_flatness, double flatness, int* out_count);
double* bz_polycurve_polyline_params(const BzPolyCurve* p, int use_flatness, double flatness, int* out_count);

double bz_polycurve_length(const BzPolyCurve* p);
double bz_polycurve_length_to(const BzPolyCurve* p, double t);
double bz_polycurve_length_between(const BzPolyCurve* p, double t1, double t2);
double bz_polycurve_step(const BzPolyCurve* p, double t, double ds);

void bz_polycurve_end_points(const BzPolyCurve* p, double out_first[2], double out_second[2]);
double* bz_polycurve_control_points(const BzPolyCurve* p, int* out_count);
void bz_polycurve_set_control_point(BzPolyCurve* p, unsigned idx, const double point[2]);

void bz_polycurve_value_at(const BzPolyCurve* p, double t, double out[2]);
double* bz_polycurve_value_at_many(const BzPolyCurve* p, const double* t, int n, int* out_count);
double bz_polycurve_curvature_at(const BzPolyCurve* p, double t);
double bz_polycurve_curvature_derivative_at(const BzPolyCurve* p, double t);
void bz_polycurve_tangent_at(const BzPolyCurve* p, double t, double out[2]);
void bz_polycurve_normal_at(const BzPolyCurve* p, double t, double out[2]);
void bz_polycurve_derivative_at(const BzPolyCurve* p, unsigned n, double t, double out[2]); /* n>=1 */
void bz_polycurve_bounding_box(const BzPolyCurve* p, double out_min[2], double out_max[2]);

double* bz_polycurve_intersections_curve(const BzPolyCurve* p, const BzCurve* c, int* out_count);
double* bz_polycurve_intersections_poly(const BzPolyCurve* p, const BzPolyCurve* other, int* out_count);
double bz_polycurve_project_point(const BzPolyCurve* p, const double point[2]);
double* bz_polycurve_project_points(const BzPolyCurve* p, const double* xy, int n_points, int* out_count);
double bz_polycurve_distance(const BzPolyCurve* p, const double point[2]);
double* bz_polycurve_distances(const BzPolyCurve* p, const double* xy, int n_points, int* out_count);

/* ============================ Utils ============================ */

/* indices into the input polyline; returned as count ints, free with bz_free */
int* bz_utils_visvalingam_wyatt(const double* xy, int n_points, int* out_count);
double* bz_utils_solve_polynomial(const double* coeffs, int n, int* out_count);
BzCurve* bz_utils_fit_bezier(const double* xy, int n_points, unsigned order);

#ifdef __cplusplus
}
#endif

#endif /* BEZIER_C_H */
