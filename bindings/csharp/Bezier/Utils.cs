// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using System;

namespace Bezier;

/// <summary>Free functions from Bezier::Utils.</summary>
public static class Utils
{
    /// <summary>Visvalingam-Whyatt simplification; returns indices into the input polyline.</summary>
    public static int[] VisvalingamWyatt(Point[] polyline)
    {
        IntPtr p = Native.bz_utils_visvalingam_wyatt(Native.Flatten(polyline), polyline.Length, out int n);
        Native.Check();
        return Native.ReadInts(p, n);
    }

    /// <summary>Real roots of a polynomial given by its coefficients.</summary>
    public static double[] SolvePolynomial(double[] coefficients)
    {
        IntPtr p = Native.bz_utils_solve_polynomial(coefficients, coefficients.Length, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    /// <summary>Fit a Bezier curve of the given order to a set of points.</summary>
    public static Curve FitBezier(Point[] points, uint order)
    {
        IntPtr ptr = Native.bz_utils_fit_bezier(Native.Flatten(points), points.Length, order);
        Native.Check();
        return new Curve(new CurveHandle(ptr));
    }
}
