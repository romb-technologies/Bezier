// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
//
// WASM module loader + heap marshalling helpers (the JS analog of the C# Native.cs).

// @ts-ignore — emitted by Emscripten at build time, resolved from dist/ at runtime.
import createBezierModule from "./bezier_wasm.mjs";

/** A point (or vector) in the xy plane. */
export interface Point {
  x: number;
  y: number;
}

/** Raised when the underlying C++ library reports an error (e.g. parameter out of range). */
export class BezierError extends Error {}

export interface BezierModule {
  _malloc(size: number): number;
  _free(ptr: number): void;
  HEAPF64: Float64Array;
  HEAP32: Int32Array;
  UTF8ToString(ptr: number): string;
  // bz_* shim functions (all exported via -sEXPORT_ALL)
  [fn: string]: any;
}

let M: BezierModule | null = null;

/** Load the WASM module once. Must be awaited before constructing curves. */
export async function load(): Promise<BezierModule> {
  if (!M) {
    const factory = createBezierModule as (opts?: Record<string, unknown>) => Promise<BezierModule>;
    M = await factory();
  }
  return M;
}

export function mod(): BezierModule {
  if (!M) throw new BezierError("Bezier WASM not loaded — call loadBezier() first");
  return M;
}

export function checkError(): void {
  const m = mod();
  const ptr = m._bz_last_error();
  if (ptr !== 0) throw new BezierError(m.UTF8ToString(ptr));
}

export function flatten(points: Point[]): Float64Array {
  const a = new Float64Array(points.length * 2);
  for (let k = 0; k < points.length; k++) {
    a[2 * k] = points[k].x;
    a[2 * k + 1] = points[k].y;
  }
  return a;
}

/** Copy a double array into the heap, run fn with (ptr, length), then free. */
export function withDoubles<T>(arr: ArrayLike<number>, fn: (ptr: number, n: number) => T): T {
  const m = mod();
  const n = arr.length;
  const ptr = m._malloc(Math.max(n, 1) * 8);
  m.HEAPF64.set(arr as unknown as ArrayLike<number>, ptr / 8);
  try {
    return fn(ptr, n);
  } finally {
    m._free(ptr);
  }
}

/** Read a single (x, y) point written by the native call into a heap buffer. */
export function readPoint(call: (outPtr: number) => void): Point {
  const m = mod();
  const ptr = m._malloc(16);
  try {
    call(ptr);
    checkError();
    return { x: m.HEAPF64[ptr / 8], y: m.HEAPF64[ptr / 8 + 1] };
  } finally {
    m._free(ptr);
  }
}

// --- variable-length returns: native returns a malloc'd pointer + writes count via int* ---

function withCount<T>(call: (countPtr: number) => number, read: (dataPtr: number, n: number) => T): T {
  const m = mod();
  const countPtr = m._malloc(4);
  try {
    const dataPtr = call(countPtr);
    checkError();
    const n = m.HEAP32[countPtr / 4];
    const result = read(dataPtr, n > 0 ? n : 0);
    if (dataPtr !== 0) m._bz_free(dataPtr);
    return result;
  } finally {
    m._free(countPtr);
  }
}

export function readDoubles(call: (countPtr: number) => number): number[] {
  return withCount(call, (ptr, n) => (n === 0 ? [] : Array.from(mod().HEAPF64.subarray(ptr / 8, ptr / 8 + n))));
}

export function readPoints(call: (countPtr: number) => number): Point[] {
  return withCount(call, (ptr, n) => {
    const pts: Point[] = [];
    const h = mod().HEAPF64;
    for (let k = 0; k < n; k++) pts.push({ x: h[ptr / 8 + 2 * k], y: h[ptr / 8 + 2 * k + 1] });
    return pts;
  });
}

export function readInts(call: (countPtr: number) => number): number[] {
  return withCount(call, (ptr, n) => (n === 0 ? [] : Array.from(mod().HEAP32.subarray(ptr / 4, ptr / 4 + n))));
}

/** Read an array of handle pointers (wasm32: 4 bytes each). */
export function readHandles(call: (countPtr: number) => number): number[] {
  return withCount(call, (ptr, n) => (n === 0 ? [] : Array.from(mod().HEAP32.subarray(ptr / 4, ptr / 4 + n))));
}
