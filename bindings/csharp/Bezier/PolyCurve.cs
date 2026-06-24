// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using System;

namespace Bezier;

/// <summary>
/// A chain of Bezier curves joined with at least C0 continuity. Wraps the native Bezier::PolyCurve.
/// The parameter t addresses subcurve n when t is in [n-1, n).
/// </summary>
public sealed class PolyCurve : IDisposable
{
    internal readonly PolyCurveHandle Handle;

    internal PolyCurve(PolyCurveHandle handle) => Handle = handle;

    /// <summary>Create an empty polycurve.</summary>
    public PolyCurve()
    {
        IntPtr ptr = Native.bz_polycurve_new();
        Native.Check();
        Handle = new PolyCurveHandle(ptr);
    }

    /// <summary>Create a polycurve from a sequence of curves (each is copied).</summary>
    public PolyCurve(params Curve[] curves)
    {
        var raw = new IntPtr[curves.Length];
        for (int k = 0; k < curves.Length; k++)
            raw[k] = curves[k].Handle.DangerousGetHandle();
        IntPtr ptr = Native.bz_polycurve_new_from(raw, curves.Length);
        Native.Check();
        GC.KeepAlive(curves);
        Handle = new PolyCurveHandle(ptr);
    }

    public PolyCurve Clone()
    {
        IntPtr p = Native.bz_polycurve_copy(Handle);
        Native.Check();
        return new PolyCurve(new PolyCurveHandle(p));
    }

    public void InsertAt(uint idx, Curve curve)
    {
        Native.bz_polycurve_insert_at(Handle, idx, curve.Handle);
        Native.Check();
    }

    public void InsertFront(Curve curve)
    {
        Native.bz_polycurve_insert_front(Handle, curve.Handle);
        Native.Check();
    }

    public void InsertBack(Curve curve)
    {
        Native.bz_polycurve_insert_back(Handle, curve.Handle);
        Native.Check();
    }

    public void RemoveAt(uint idx)
    {
        Native.bz_polycurve_remove_at(Handle, idx);
        Native.Check();
    }

    public void RemoveFirst()
    {
        Native.bz_polycurve_remove_first(Handle);
        Native.Check();
    }

    public void RemoveBack()
    {
        Native.bz_polycurve_remove_back(Handle);
        Native.Check();
    }

    public uint Size()
    {
        uint v = Native.bz_polycurve_size(Handle);
        Native.Check();
        return v;
    }

    public uint CurveIdx(double t)
    {
        uint v = Native.bz_polycurve_curve_idx(Handle, t);
        Native.Check();
        return v;
    }

    /// <summary>The subcurve at idx, as an independent copy.</summary>
    public Curve GetCurve(uint idx)
    {
        IntPtr p = Native.bz_polycurve_curve(Handle, idx);
        Native.Check();
        return new Curve(new CurveHandle(p));
    }

    public Point[] Polyline()
    {
        IntPtr p = Native.bz_polycurve_polyline(Handle, 0, 0.0, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public Point[] Polyline(double flatness)
    {
        IntPtr p = Native.bz_polycurve_polyline(Handle, 1, flatness, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public double[] PolylineParams()
    {
        IntPtr p = Native.bz_polycurve_polyline_params(Handle, 0, 0.0, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public double[] PolylineParams(double flatness)
    {
        IntPtr p = Native.bz_polycurve_polyline_params(Handle, 1, flatness, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public double Length()
    {
        double v = Native.bz_polycurve_length(Handle);
        Native.Check();
        return v;
    }

    public double Length(double t)
    {
        double v = Native.bz_polycurve_length_to(Handle, t);
        Native.Check();
        return v;
    }

    public double Length(double t1, double t2)
    {
        double v = Native.bz_polycurve_length_between(Handle, t1, t2);
        Native.Check();
        return v;
    }

    public double Step(double t, double ds)
    {
        double v = Native.bz_polycurve_step(Handle, t, ds);
        Native.Check();
        return v;
    }

    public (Point First, Point Last) EndPoints()
    {
        var a = new double[2];
        var b = new double[2];
        Native.bz_polycurve_end_points(Handle, a, b);
        Native.Check();
        return (new Point(a[0], a[1]), new Point(b[0], b[1]));
    }

    public Point[] ControlPoints()
    {
        IntPtr p = Native.bz_polycurve_control_points(Handle, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public void SetControlPoint(uint idx, Point point)
    {
        Native.bz_polycurve_set_control_point(Handle, idx, new[] { point.X, point.Y });
        Native.Check();
    }

    public Point ValueAt(double t)
    {
        var b = new double[2];
        Native.bz_polycurve_value_at(Handle, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public Point[] ValueAt(double[] t)
    {
        IntPtr p = Native.bz_polycurve_value_at_many(Handle, t, t.Length, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public double CurvatureAt(double t)
    {
        double v = Native.bz_polycurve_curvature_at(Handle, t);
        Native.Check();
        return v;
    }

    public double CurvatureDerivativeAt(double t)
    {
        double v = Native.bz_polycurve_curvature_derivative_at(Handle, t);
        Native.Check();
        return v;
    }

    public Point TangentAt(double t)
    {
        var b = new double[2];
        Native.bz_polycurve_tangent_at(Handle, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public Point NormalAt(double t)
    {
        var b = new double[2];
        Native.bz_polycurve_normal_at(Handle, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public Point DerivativeAt(double t, uint n = 1)
    {
        var b = new double[2];
        Native.bz_polycurve_derivative_at(Handle, n, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public (Point Min, Point Max) BoundingBox()
    {
        var lo = new double[2];
        var hi = new double[2];
        Native.bz_polycurve_bounding_box(Handle, lo, hi);
        Native.Check();
        return (new Point(lo[0], lo[1]), new Point(hi[0], hi[1]));
    }

    public Point[] Intersections(Curve curve)
    {
        IntPtr p = Native.bz_polycurve_intersections_curve(Handle, curve.Handle, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public Point[] Intersections(PolyCurve other)
    {
        IntPtr p = Native.bz_polycurve_intersections_poly(Handle, other.Handle, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public double ProjectPoint(Point point)
    {
        double v = Native.bz_polycurve_project_point(Handle, new[] { point.X, point.Y });
        Native.Check();
        return v;
    }

    public double[] ProjectPoint(Point[] points)
    {
        IntPtr p = Native.bz_polycurve_project_points(Handle, Native.Flatten(points), points.Length, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public double Distance(Point point)
    {
        double v = Native.bz_polycurve_distance(Handle, new[] { point.X, point.Y });
        Native.Check();
        return v;
    }

    public double[] Distance(Point[] points)
    {
        IntPtr p = Native.bz_polycurve_distances(Handle, Native.Flatten(points), points.Length, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public override string ToString() => $"<Bezier.PolyCurve size={Size()}>";

    public void Dispose() => Handle.Dispose();
}
