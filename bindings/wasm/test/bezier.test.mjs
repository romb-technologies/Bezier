// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
import test from "node:test";
import assert from "node:assert/strict";
import { loadBezier, BezierError } from "../dist/index.js";

const B = await loadBezier();
const cp = () => [
  { x: 0, y: 0 },
  { x: 1, y: 2 },
  { x: 3, y: 3 },
  { x: 4, y: 0 },
];

const close = (a, b, tol = 1e-9) => Math.abs(a.x - b.x) < tol && Math.abs(a.y - b.y) < tol;

test("construct and order", () => {
  const c = new B.Curve(cp());
  assert.equal(c.order(), 3);
  assert.equal(c.controlPoints().length, 4);
  c.dispose();
});

test("endpoints match valueAt", () => {
  const c = new B.Curve(cp());
  assert.ok(close(c.valueAt(0), cp()[0]));
  assert.ok(close(c.valueAt(1), cp()[3]));
  c.dispose();
});

test("batch valueAt count", () => {
  const c = new B.Curve(cp());
  const pts = c.valueAt([0, 0.5, 1]);
  assert.equal(pts.length, 3);
  c.dispose();
});

test("length positive", () => {
  const c = new B.Curve(cp());
  assert.ok(c.length() > 0);
  c.dispose();
});

test("derivative lowers order", () => {
  const c = new B.Curve(cp());
  const d1 = c.derivative();
  const d2 = c.derivative(2);
  assert.equal(d1.order(), 2);
  assert.equal(d2.order(), 1);
  c.dispose();
  d1.dispose();
  d2.dispose();
});

test("split roundtrips endpoints", () => {
  const c = new B.Curve(cp());
  const [left, right] = c.splitCurve(0.5);
  assert.ok(close(left.valueAt(0), cp()[0]));
  assert.ok(close(right.valueAt(1), cp()[3]));
  assert.ok(close(left.valueAt(1), right.valueAt(0))); // C0 at the seam
  c.dispose();
  left.dispose();
  right.dispose();
});

test("bounding box ordered", () => {
  const c = new B.Curve(cp());
  const { min, max } = c.boundingBox();
  assert.ok(min.x <= max.x && min.y <= max.y);
  c.dispose();
});

test("polycurve size", () => {
  const c1 = new B.Curve(cp());
  const c2 = new B.Curve(cp().map((p) => ({ x: p.x + 4, y: p.y + 4 })));
  const pc = new B.PolyCurve(c1, c2);
  assert.equal(pc.size(), 2);
  const sub = pc.getCurve(0);
  assert.equal(sub.order(), 3);
  c1.dispose();
  c2.dispose();
  sub.dispose();
  pc.dispose();
});

test("out-of-range param throws, not NaN", () => {
  const c = new B.Curve(cp());
  assert.throws(() => c.length(5.0), BezierError);
  c.dispose();
});
