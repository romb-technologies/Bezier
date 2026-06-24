# Bezier — WebAssembly (TS/JS) bindings

WebAssembly bindings for the [Bezier](../../README.md) C++ library. The shared C-ABI shim
(`../c-abi/`) and the core sources are compiled to WASM with [Emscripten](https://emscripten.org),
and wrapped in a small TypeScript layer. Ships as an ESM npm package.

## Build

Needs [emsdk](https://emscripten.org/docs/getting_started/downloads.html) (Emscripten),
Node.js, CMake, and Eigen3. Use a **recent Emscripten** (tested with 6.0.1) — older versions
(e.g. 3.1.6) are too old: their libc++ rejects over-aligned `std::vector<Eigen::Vector2d>`
and the ESM output doesn't load under Node. Install emsdk and `source emsdk_env.sh`.

```bash
# 1. compile the C++ to WASM -> dist/bezier_wasm.{mjs,wasm}
source ~/emsdk/emsdk_env.sh
emcmake cmake -S bindings/wasm -B bindings/wasm/build
cmake --build bindings/wasm/build

# 2. build the TS wrapper and run the tests
cd bindings/wasm
npm install
npm run build      # tsc -> dist/index.js + .d.ts
npm test           # node --test
```

## Usage

```ts
import { loadBezier } from "@romb/bezier";

const { Curve, PolyCurve } = await loadBezier();   // WASM loads asynchronously

const c = new Curve([
  { x: 0, y: 0 }, { x: 1, y: 2 }, { x: 3, y: 3 }, { x: 4, y: 0 },
]);

console.log(c.order());           // 3
console.log(c.valueAt(0.5));      // { x, y }
console.log(c.valueAt([0, 0.5, 1])); // [{x,y}, ...]
console.log(c.length());
const d = c.derivative();         // order-2 curve

const pc = new PolyCurve(c, Curve.fromPolyline(c.controlPoints()));
console.log(pc.size());           // 2

// objects own native memory — dispose when done
d.dispose();
c.dispose();
pc.dispose();
```

Method names mirror the C++ / Python / C# bindings. Points are `{ x, y }` objects;
`boundingBox()` returns `{ min, max }`. Out-of-range parameters and other C++ errors are
thrown as `BezierError` (never NaN). The free functions `visvalingamWyatt`,
`solvePolynomial`, and `fitBezier` live on the `Utils` export.

## Memory

WebAssembly has no GC over the C++ objects, so each `Curve`/`PolyCurve` (including those
returned by `derivative()`, `splitCurve()`, `getCurve()`, and the static factories) owns
native memory. Call `dispose()` when finished. A `FinalizationRegistry` frees leaked objects
on a best-effort basis, but explicit `dispose()` is the reliable path.
