# Bezier — C# bindings

C# bindings for the [Bezier](../../README.md) C++ library, via a thin C ABI shim
(`native/`) called from idiomatic C# wrappers over P/Invoke. Targets .NET 8, linux-x64.

## Build

Needs a C++17 compiler, CMake, Eigen3, and the .NET 8 SDK (`dotnet-sdk-8.0`).

```bash
# 1. build the native shim -> libbezier_c.so
cmake -S bindings/csharp/native -B bindings/csharp/native/build
cmake --build bindings/csharp/native/build

# 2. build / test the managed library (copies the .so next to the assembly)
dotnet test bindings/csharp/Bezier.sln
```

## Usage

```csharp
using Bezier;

using var c = new Curve(new[]
{
    new Point(0, 0), new Point(1, 2), new Point(3, 3), new Point(4, 0),
});

Console.WriteLine(c.Order());            // 3
Console.WriteLine(c.ValueAt(0.5));       // point on the curve
Point[] pts = c.ValueAt(new[] { 0.0, 0.5, 1.0 });
Console.WriteLine(c.Length());
using Curve d = c.Derivative();          // order-2 curve (independent copy)

using var pc = new PolyCurve(c, Curve.FromPolyline(c.ControlPoints()));
Console.WriteLine(pc.Size());            // 2
```

Method names mirror the C++ API (PascalCase). Points are `Bezier.Point` structs;
`BoundingBox()` returns a `(Point Min, Point Max)` tuple. Out-of-range parameters and other
C++ errors are surfaced as `BezierException` (never NaN). `Curve`/`PolyCurve` are
`IDisposable` — wrap them in `using`. The free functions `VisvalingamWyatt`,
`SolvePolynomial`, and `FitBezier` live in the static `Bezier.Utils` class.

## Notes

- `Derivative()` and `PolyCurve.GetCurve()` return independent copies, not views into the
  parent's internal cache — safe to keep and dispose on their own.
- Native library resolution relies on `libbezier_c.so` sitting next to the managed assembly
  (handled by the `.csproj` after step 1). For deployment, package it under
  `runtimes/linux-x64/native/` in a NuGet package (not done here).
