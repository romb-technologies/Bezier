// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
//
// Public TS API for the Bezier WASM bindings. Method names mirror the C++ / Python / C#
// bindings. Objects own native memory: call dispose() when done (a FinalizationRegistry is a
// best-effort safety net, but WASM has no GC over C++ objects).

import {
  BezierError,
  Point,
  checkError,
  flatten,
  load,
  mod,
  readDoubles,
  readHandles,
  readInts,
  readPoint,
  readPoints,
  withDoubles,
} from "./runtime.js";

export { BezierError, Point };

const curveRegistry = new FinalizationRegistry<number>((ptr) => mod()._bz_curve_destroy(ptr));
const polyRegistry = new FinalizationRegistry<number>((ptr) => mod()._bz_polycurve_destroy(ptr));

export class Curve {
  /** @internal native pointer */
  ptr: number;

  constructor(controlPoints: Point[]) {
    this.ptr = withDoubles(flatten(controlPoints), (p) => mod()._bz_curve_new(p, controlPoints.length));
    checkError();
    curveRegistry.register(this, this.ptr, this);
  }

  /** @internal wrap an existing native curve handle */
  static _wrap(ptr: number): Curve {
    const c: Curve = Object.create(Curve.prototype);
    c.ptr = ptr;
    curveRegistry.register(c, ptr, c);
    return c;
  }

  private scalar(v: number): number {
    checkError();
    return v;
  }

  clone(): Curve {
    return Curve._wrap(this.scalarPtr(mod()._bz_curve_copy(this.ptr)));
  }

  private scalarPtr(ptr: number): number {
    checkError();
    return ptr;
  }

  order(): number {
    return this.scalar(mod()._bz_curve_order(this.ptr));
  }

  controlPoints(): Point[] {
    return readPoints((cp) => mod()._bz_curve_control_points(this.ptr, cp));
  }

  controlPoint(idx: number): Point {
    return readPoint((out) => mod()._bz_curve_control_point(this.ptr, idx, out));
  }

  endPoints(): [Point, Point] {
    const m = mod();
    const a = m._malloc(16);
    const b = m._malloc(16);
    try {
      m._bz_curve_end_points(this.ptr, a, b);
      checkError();
      return [
        { x: m.HEAPF64[a / 8], y: m.HEAPF64[a / 8 + 1] },
        { x: m.HEAPF64[b / 8], y: m.HEAPF64[b / 8 + 1] },
      ];
    } finally {
      m._free(a);
      m._free(b);
    }
  }

  polyline(flatness?: number): Point[] {
    const use = flatness === undefined ? 0 : 1;
    return readPoints((cp) => mod()._bz_curve_polyline(this.ptr, use, flatness ?? 0, cp));
  }

  polylineParams(flatness?: number): number[] {
    const use = flatness === undefined ? 0 : 1;
    return readDoubles((cp) => mod()._bz_curve_polyline_params(this.ptr, use, flatness ?? 0, cp));
  }

  length(t1?: number, t2?: number): number {
    const m = mod();
    if (t1 === undefined) return this.scalar(m._bz_curve_length(this.ptr));
    if (t2 === undefined) return this.scalar(m._bz_curve_length_to(this.ptr, t1));
    return this.scalar(m._bz_curve_length_between(this.ptr, t1, t2));
  }

  step(t: number, ds: number): number {
    return this.scalar(mod()._bz_curve_step(this.ptr, t, ds));
  }

  reverse(): void {
    mod()._bz_curve_reverse(this.ptr);
    checkError();
  }

  setControlPoint(idx: number, point: Point): void {
    withDoubles([point.x, point.y], (p) => mod()._bz_curve_set_control_point(this.ptr, idx, p));
    checkError();
  }

  raiseOrder(): void {
    mod()._bz_curve_raise_order(this.ptr);
    checkError();
  }

  lowerOrder(): void {
    mod()._bz_curve_lower_order(this.ptr);
    checkError();
  }

  valueAt(t: number): Point;
  valueAt(t: number[]): Point[];
  valueAt(t: number | number[]): Point | Point[] {
    if (Array.isArray(t)) {
      return withDoubles(t, (tp, n) => readPoints((cp) => mod()._bz_curve_value_at_many(this.ptr, tp, n, cp)));
    }
    return readPoint((out) => mod()._bz_curve_value_at(this.ptr, t, out));
  }

  curvatureAt(t: number): number {
    return this.scalar(mod()._bz_curve_curvature_at(this.ptr, t));
  }

  curvatureDerivativeAt(t: number): number {
    return this.scalar(mod()._bz_curve_curvature_derivative_at(this.ptr, t));
  }

  tangentAt(t: number): Point {
    return readPoint((out) => mod()._bz_curve_tangent_at(this.ptr, t, out));
  }

  normalAt(t: number): Point {
    return readPoint((out) => mod()._bz_curve_normal_at(this.ptr, t, out));
  }

  /** The n-th derivative as an independent curve (default first derivative). */
  derivative(n = 1): Curve {
    return Curve._wrap(this.scalarPtr(mod()._bz_curve_derivative(this.ptr, n)));
  }

  derivativeAt(t: number, n = 1): Point {
    return readPoint((out) => mod()._bz_curve_derivative_at(this.ptr, n, t, out));
  }

  roots(): number[] {
    return readDoubles((cp) => mod()._bz_curve_roots(this.ptr, cp));
  }

  extrema(): number[] {
    return readDoubles((cp) => mod()._bz_curve_extrema(this.ptr, cp));
  }

  boundingBox(): { min: Point; max: Point } {
    const m = mod();
    const lo = m._malloc(16);
    const hi = m._malloc(16);
    try {
      m._bz_curve_bounding_box(this.ptr, lo, hi);
      checkError();
      return {
        min: { x: m.HEAPF64[lo / 8], y: m.HEAPF64[lo / 8 + 1] },
        max: { x: m.HEAPF64[hi / 8], y: m.HEAPF64[hi / 8 + 1] },
      };
    } finally {
      m._free(lo);
      m._free(hi);
    }
  }

  splitCurve(t?: number): [Curve, Curve];
  splitCurve(t: number[]): Curve[];
  splitCurve(t: number | number[] = 0.5): [Curve, Curve] | Curve[] {
    const m = mod();
    if (Array.isArray(t)) {
      const handles = withDoubles(t, (tp, n) => readHandles((cp) => m._bz_curve_split(this.ptr, tp, n, cp)));
      return handles.map((h) => Curve._wrap(h));
    }
    const left = m._malloc(4);
    const right = m._malloc(4);
    try {
      m._bz_curve_split_at(this.ptr, t, left, right);
      checkError();
      return [Curve._wrap(m.HEAP32[left / 4]), Curve._wrap(m.HEAP32[right / 4])];
    } finally {
      m._free(left);
      m._free(right);
    }
  }

  intersections(other: Curve): Point[] {
    return readPoints((cp) => mod()._bz_curve_intersections(this.ptr, other.ptr, cp));
  }

  projectPoint(point: Point): number {
    const v = withDoubles([point.x, point.y], (p) => mod()._bz_curve_project_point(this.ptr, p));
    checkError();
    return v;
  }

  distance(point: Point): number {
    const v = withDoubles([point.x, point.y], (p) => mod()._bz_curve_distance(this.ptr, p));
    checkError();
    return v;
  }

  applyContinuity(other: Curve, betaCoeffs: number[]): void {
    withDoubles(betaCoeffs, (b, n) => mod()._bz_curve_apply_continuity(this.ptr, other.ptr, b, n));
    checkError();
  }

  static offsetCurve(curve: Curve, offset: number, order = 0): Curve {
    const ptr = mod()._bz_curve_offset(curve.ptr, offset, order);
    checkError();
    return Curve._wrap(ptr);
  }

  static joinCurves(curve1: Curve, curve2: Curve, order = 0): Curve {
    const ptr = mod()._bz_curve_join(curve1.ptr, curve2.ptr, order);
    checkError();
    return Curve._wrap(ptr);
  }

  static fromPolyline(polyline: Point[], order = 0): Curve {
    const ptr = withDoubles(flatten(polyline), (p) => mod()._bz_curve_from_polyline(p, polyline.length, order));
    checkError();
    return Curve._wrap(ptr);
  }

  dispose(): void {
    if (this.ptr !== 0) {
      curveRegistry.unregister(this);
      mod()._bz_curve_destroy(this.ptr);
      this.ptr = 0;
    }
  }

  toString(): string {
    return `<Bezier.Curve order=${this.order()}>`;
  }
}

export class PolyCurve {
  /** @internal native pointer */
  ptr: number;

  constructor(...curves: Curve[]) {
    const m = mod();
    if (curves.length === 0) {
      this.ptr = m._bz_polycurve_new();
      checkError();
    } else {
      const arr = m._malloc(curves.length * 4);
      try {
        for (let k = 0; k < curves.length; k++) m.HEAP32[arr / 4 + k] = curves[k].ptr;
        this.ptr = m._bz_polycurve_new_from(arr, curves.length);
        checkError();
      } finally {
        m._free(arr);
      }
    }
    polyRegistry.register(this, this.ptr, this);
  }

  /** @internal */
  static _wrap(ptr: number): PolyCurve {
    const p: PolyCurve = Object.create(PolyCurve.prototype);
    p.ptr = ptr;
    polyRegistry.register(p, ptr, p);
    return p;
  }

  private scalar(v: number): number {
    checkError();
    return v;
  }

  insertAt(idx: number, curve: Curve): void {
    mod()._bz_polycurve_insert_at(this.ptr, idx, curve.ptr);
    checkError();
  }

  insertFront(curve: Curve): void {
    mod()._bz_polycurve_insert_front(this.ptr, curve.ptr);
    checkError();
  }

  insertBack(curve: Curve): void {
    mod()._bz_polycurve_insert_back(this.ptr, curve.ptr);
    checkError();
  }

  removeAt(idx: number): void {
    mod()._bz_polycurve_remove_at(this.ptr, idx);
    checkError();
  }

  removeFirst(): void {
    mod()._bz_polycurve_remove_first(this.ptr);
    checkError();
  }

  removeBack(): void {
    mod()._bz_polycurve_remove_back(this.ptr);
    checkError();
  }

  size(): number {
    return this.scalar(mod()._bz_polycurve_size(this.ptr));
  }

  curveIdx(t: number): number {
    return this.scalar(mod()._bz_polycurve_curve_idx(this.ptr, t));
  }

  /** The subcurve at idx, as an independent copy. */
  getCurve(idx: number): Curve {
    const ptr = mod()._bz_polycurve_curve(this.ptr, idx);
    checkError();
    return Curve._wrap(ptr);
  }

  polyline(flatness?: number): Point[] {
    const use = flatness === undefined ? 0 : 1;
    return readPoints((cp) => mod()._bz_polycurve_polyline(this.ptr, use, flatness ?? 0, cp));
  }

  polylineParams(flatness?: number): number[] {
    const use = flatness === undefined ? 0 : 1;
    return readDoubles((cp) => mod()._bz_polycurve_polyline_params(this.ptr, use, flatness ?? 0, cp));
  }

  length(t1?: number, t2?: number): number {
    const m = mod();
    if (t1 === undefined) return this.scalar(m._bz_polycurve_length(this.ptr));
    if (t2 === undefined) return this.scalar(m._bz_polycurve_length_to(this.ptr, t1));
    return this.scalar(m._bz_polycurve_length_between(this.ptr, t1, t2));
  }

  step(t: number, ds: number): number {
    return this.scalar(mod()._bz_polycurve_step(this.ptr, t, ds));
  }

  endPoints(): [Point, Point] {
    const m = mod();
    const a = m._malloc(16);
    const b = m._malloc(16);
    try {
      m._bz_polycurve_end_points(this.ptr, a, b);
      checkError();
      return [
        { x: m.HEAPF64[a / 8], y: m.HEAPF64[a / 8 + 1] },
        { x: m.HEAPF64[b / 8], y: m.HEAPF64[b / 8 + 1] },
      ];
    } finally {
      m._free(a);
      m._free(b);
    }
  }

  controlPoints(): Point[] {
    return readPoints((cp) => mod()._bz_polycurve_control_points(this.ptr, cp));
  }

  setControlPoint(idx: number, point: Point): void {
    withDoubles([point.x, point.y], (p) => mod()._bz_polycurve_set_control_point(this.ptr, idx, p));
    checkError();
  }

  valueAt(t: number): Point;
  valueAt(t: number[]): Point[];
  valueAt(t: number | number[]): Point | Point[] {
    if (Array.isArray(t)) {
      return withDoubles(t, (tp, n) => readPoints((cp) => mod()._bz_polycurve_value_at_many(this.ptr, tp, n, cp)));
    }
    return readPoint((out) => mod()._bz_polycurve_value_at(this.ptr, t, out));
  }

  curvatureAt(t: number): number {
    return this.scalar(mod()._bz_polycurve_curvature_at(this.ptr, t));
  }

  curvatureDerivativeAt(t: number): number {
    return this.scalar(mod()._bz_polycurve_curvature_derivative_at(this.ptr, t));
  }

  tangentAt(t: number): Point {
    return readPoint((out) => mod()._bz_polycurve_tangent_at(this.ptr, t, out));
  }

  normalAt(t: number): Point {
    return readPoint((out) => mod()._bz_polycurve_normal_at(this.ptr, t, out));
  }

  derivativeAt(t: number, n = 1): Point {
    return readPoint((out) => mod()._bz_polycurve_derivative_at(this.ptr, n, t, out));
  }

  boundingBox(): { min: Point; max: Point } {
    const m = mod();
    const lo = m._malloc(16);
    const hi = m._malloc(16);
    try {
      m._bz_polycurve_bounding_box(this.ptr, lo, hi);
      checkError();
      return {
        min: { x: m.HEAPF64[lo / 8], y: m.HEAPF64[lo / 8 + 1] },
        max: { x: m.HEAPF64[hi / 8], y: m.HEAPF64[hi / 8 + 1] },
      };
    } finally {
      m._free(lo);
      m._free(hi);
    }
  }

  intersections(target: Curve | PolyCurve): Point[] {
    const m = mod();
    return readPoints((cp) =>
      target instanceof Curve
        ? m._bz_polycurve_intersections_curve(this.ptr, target.ptr, cp)
        : m._bz_polycurve_intersections_poly(this.ptr, target.ptr, cp),
    );
  }

  projectPoint(point: Point): number;
  projectPoint(points: Point[]): number[];
  projectPoint(point: Point | Point[]): number | number[] {
    const m = mod();
    if (Array.isArray(point)) {
      return withDoubles(flatten(point), (p) =>
        readDoubles((cp) => m._bz_polycurve_project_points(this.ptr, p, point.length, cp)),
      );
    }
    const v = withDoubles([point.x, point.y], (p) => m._bz_polycurve_project_point(this.ptr, p));
    checkError();
    return v;
  }

  distance(point: Point): number;
  distance(points: Point[]): number[];
  distance(point: Point | Point[]): number | number[] {
    const m = mod();
    if (Array.isArray(point)) {
      return withDoubles(flatten(point), (p) =>
        readDoubles((cp) => m._bz_polycurve_distances(this.ptr, p, point.length, cp)),
      );
    }
    const v = withDoubles([point.x, point.y], (p) => m._bz_polycurve_distance(this.ptr, p));
    checkError();
    return v;
  }

  dispose(): void {
    if (this.ptr !== 0) {
      polyRegistry.unregister(this);
      mod()._bz_polycurve_destroy(this.ptr);
      this.ptr = 0;
    }
  }

  toString(): string {
    return `<Bezier.PolyCurve size=${this.size()}>`;
  }
}

export const Utils = {
  visvalingamWyatt(polyline: Point[]): number[] {
    // returns indices into the input polyline (ints)
    return withDoubles(flatten(polyline), (p) =>
      readInts((cp) => mod()._bz_utils_visvalingam_wyatt(p, polyline.length, cp)),
    );
  },

  solvePolynomial(coefficients: number[]): number[] {
    return withDoubles(coefficients, (p, n) =>
      readDoubles((cp) => mod()._bz_utils_solve_polynomial(p, n, cp)),
    );
  },

  fitBezier(points: Point[], order: number): Curve {
    const ptr = withDoubles(flatten(points), (p) => mod()._bz_utils_fit_bezier(p, points.length, order));
    checkError();
    return Curve._wrap(ptr);
  },
};

/** Load the WASM module and return the API bound to it. Call once before use. */
export async function loadBezier(): Promise<{
  Curve: typeof Curve;
  PolyCurve: typeof PolyCurve;
  Utils: typeof Utils;
}> {
  await load();
  return { Curve, PolyCurve, Utils };
}
