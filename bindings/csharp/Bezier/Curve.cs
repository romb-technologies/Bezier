// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using System;

namespace Bezier;

/// <summary>A Bezier curve of any order. Wraps the native Bezier::Curve.</summary>
public sealed class Curve : IDisposable
{
    internal readonly CurveHandle Handle;

    internal Curve(CurveHandle handle) => Handle = handle;

    private static Curve Wrap(IntPtr ptr)
    {
        Native.Check();
        return new Curve(new CurveHandle(ptr));
    }

    /// <summary>Create a curve from control points (each row a point).</summary>
    public Curve(Point[] controlPoints)
    {
        IntPtr ptr = Native.bz_curve_new(Native.Flatten(controlPoints), controlPoints.Length);
        Native.Check();
        Handle = new CurveHandle(ptr);
    }

    public Curve Clone() => Wrap(Native.bz_curve_copy(Handle));

    public uint Order()
    {
        uint v = Native.bz_curve_order(Handle);
        Native.Check();
        return v;
    }

    public Point[] ControlPoints()
    {
        IntPtr p = Native.bz_curve_control_points(Handle, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public Point ControlPoint(uint idx)
    {
        var b = new double[2];
        Native.bz_curve_control_point(Handle, idx, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public (Point First, Point Last) EndPoints()
    {
        var a = new double[2];
        var b = new double[2];
        Native.bz_curve_end_points(Handle, a, b);
        Native.Check();
        return (new Point(a[0], a[1]), new Point(b[0], b[1]));
    }

    public Point[] Polyline()
    {
        IntPtr p = Native.bz_curve_polyline(Handle, 0, 0.0, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public Point[] Polyline(double flatness)
    {
        IntPtr p = Native.bz_curve_polyline(Handle, 1, flatness, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public double[] PolylineParams()
    {
        IntPtr p = Native.bz_curve_polyline_params(Handle, 0, 0.0, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public double[] PolylineParams(double flatness)
    {
        IntPtr p = Native.bz_curve_polyline_params(Handle, 1, flatness, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public double Length()
    {
        double v = Native.bz_curve_length(Handle);
        Native.Check();
        return v;
    }

    public double Length(double t)
    {
        double v = Native.bz_curve_length_to(Handle, t);
        Native.Check();
        return v;
    }

    public double Length(double t1, double t2)
    {
        double v = Native.bz_curve_length_between(Handle, t1, t2);
        Native.Check();
        return v;
    }

    public double Step(double t, double ds)
    {
        double v = Native.bz_curve_step(Handle, t, ds);
        Native.Check();
        return v;
    }

    public void Reverse()
    {
        Native.bz_curve_reverse(Handle);
        Native.Check();
    }

    public void SetControlPoint(uint idx, Point point)
    {
        Native.bz_curve_set_control_point(Handle, idx, new[] { point.X, point.Y });
        Native.Check();
    }

    public void RaiseOrder()
    {
        Native.bz_curve_raise_order(Handle);
        Native.Check();
    }

    public void LowerOrder()
    {
        Native.bz_curve_lower_order(Handle);
        Native.Check();
    }

    public Point ValueAt(double t)
    {
        var b = new double[2];
        Native.bz_curve_value_at(Handle, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public Point[] ValueAt(double[] t)
    {
        IntPtr p = Native.bz_curve_value_at_many(Handle, t, t.Length, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public double CurvatureAt(double t)
    {
        double v = Native.bz_curve_curvature_at(Handle, t);
        Native.Check();
        return v;
    }

    public double CurvatureDerivativeAt(double t)
    {
        double v = Native.bz_curve_curvature_derivative_at(Handle, t);
        Native.Check();
        return v;
    }

    public Point TangentAt(double t)
    {
        var b = new double[2];
        Native.bz_curve_tangent_at(Handle, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public Point NormalAt(double t)
    {
        var b = new double[2];
        Native.bz_curve_normal_at(Handle, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    /// <summary>The n-th derivative as an independent curve (default first derivative).</summary>
    public Curve Derivative(uint n = 1) => Wrap(Native.bz_curve_derivative(Handle, n));

    public Point DerivativeAt(double t, uint n = 1)
    {
        var b = new double[2];
        Native.bz_curve_derivative_at(Handle, n, t, b);
        Native.Check();
        return new Point(b[0], b[1]);
    }

    public double[] Roots()
    {
        IntPtr p = Native.bz_curve_roots(Handle, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public double[] Extrema()
    {
        IntPtr p = Native.bz_curve_extrema(Handle, out int n);
        Native.Check();
        return Native.ReadDoubles(p, n);
    }

    public (Point Min, Point Max) BoundingBox()
    {
        var lo = new double[2];
        var hi = new double[2];
        Native.bz_curve_bounding_box(Handle, lo, hi);
        Native.Check();
        return (new Point(lo[0], lo[1]), new Point(hi[0], hi[1]));
    }

    public Curve[] SplitCurve(double[] t)
    {
        IntPtr p = Native.bz_curve_split(Handle, t, t.Length, out int n);
        Native.Check();
        IntPtr[] handles = Native.ReadHandles(p, n);
        var parts = new Curve[handles.Length];
        for (int k = 0; k < handles.Length; k++)
            parts[k] = new Curve(new CurveHandle(handles[k]));
        return parts;
    }

    public (Curve Left, Curve Right) SplitCurve(double t = 0.5)
    {
        Native.bz_curve_split_at(Handle, t, out IntPtr left, out IntPtr right);
        Native.Check();
        return (new Curve(new CurveHandle(left)), new Curve(new CurveHandle(right)));
    }

    public Point[] Intersections(Curve other)
    {
        IntPtr p = Native.bz_curve_intersections(Handle, other.Handle, out int n);
        Native.Check();
        return Native.ReadPoints(p, n);
    }

    public double ProjectPoint(Point point)
    {
        double v = Native.bz_curve_project_point(Handle, new[] { point.X, point.Y });
        Native.Check();
        return v;
    }

    public double Distance(Point point)
    {
        double v = Native.bz_curve_distance(Handle, new[] { point.X, point.Y });
        Native.Check();
        return v;
    }

    public void ApplyContinuity(Curve other, double[] betaCoeffs)
    {
        Native.bz_curve_apply_continuity(Handle, other.Handle, betaCoeffs, betaCoeffs.Length);
        Native.Check();
    }

    public static Curve OffsetCurve(Curve curve, double offset, uint order = 0) =>
        Wrap(Native.bz_curve_offset(curve.Handle, offset, order));

    public static Curve JoinCurves(Curve curve1, Curve curve2, uint order = 0) =>
        Wrap(Native.bz_curve_join(curve1.Handle, curve2.Handle, order));

    public static Curve FromPolyline(Point[] polyline, uint order = 0) =>
        Wrap(Native.bz_curve_from_polyline(Native.Flatten(polyline), polyline.Length, order));

    public override string ToString() => $"<Bezier.Curve order={Order()}>";

    public void Dispose() => Handle.Dispose();
}
