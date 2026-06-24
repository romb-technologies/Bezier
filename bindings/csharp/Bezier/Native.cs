// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using System;
using System.Runtime.InteropServices;

namespace Bezier;

/// <summary>P/Invoke declarations and marshalling helpers for the libbezier_c shim.</summary>
internal static partial class Native
{
    private const string Lib = "bezier_c";

    // --- error handling / memory ---
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    private static extern IntPtr bz_last_error();

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_free(IntPtr p);

    /// <summary>Throw if the most recent native call recorded an error.</summary>
    internal static void Check()
    {
        IntPtr e = bz_last_error();
        if (e != IntPtr.Zero)
            throw new BezierException(Marshal.PtrToStringUTF8(e) ?? "Bezier error");
    }

    // --- marshalling helpers (always call Check() before reading native buffers) ---
    internal static double[] Flatten(Point[] pts)
    {
        var d = new double[pts.Length * 2];
        for (int k = 0; k < pts.Length; k++)
        {
            d[2 * k] = pts[k].X;
            d[2 * k + 1] = pts[k].Y;
        }
        return d;
    }

    internal static Point[] ReadPoints(IntPtr ptr, int count)
    {
        if (count <= 0)
            return Array.Empty<Point>();
        var flat = new double[count * 2];
        Marshal.Copy(ptr, flat, 0, count * 2);
        bz_free(ptr);
        var pts = new Point[count];
        for (int k = 0; k < count; k++)
            pts[k] = new Point(flat[2 * k], flat[2 * k + 1]);
        return pts;
    }

    internal static double[] ReadDoubles(IntPtr ptr, int count)
    {
        if (count <= 0)
            return Array.Empty<double>();
        var d = new double[count];
        Marshal.Copy(ptr, d, 0, count);
        bz_free(ptr);
        return d;
    }

    internal static int[] ReadInts(IntPtr ptr, int count)
    {
        if (count <= 0)
            return Array.Empty<int>();
        var d = new int[count];
        Marshal.Copy(ptr, d, 0, count);
        bz_free(ptr);
        return d;
    }

    internal static IntPtr[] ReadHandles(IntPtr ptr, int count)
    {
        if (count <= 0)
            return Array.Empty<IntPtr>();
        var h = new IntPtr[count];
        Marshal.Copy(ptr, h, 0, count);
        bz_free(ptr);
        return h;
    }

    // ============================ Curve ============================
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_new(double[] xy, int nPoints);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_copy(CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_destroy(IntPtr c);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern uint bz_curve_order(CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_control_points(CurveHandle c, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_control_point(CurveHandle c, uint idx, [Out] double[] outPoint);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_end_points(CurveHandle c, [Out] double[] first, [Out] double[] second);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_polyline(CurveHandle c, int useFlatness, double flatness, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_polyline_params(CurveHandle c, int useFlatness, double flatness, out int count);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_length(CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_length_to(CurveHandle c, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_length_between(CurveHandle c, double t1, double t2);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_step(CurveHandle c, double t, double ds);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_reverse(CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_set_control_point(CurveHandle c, uint idx, double[] point);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_raise_order(CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_lower_order(CurveHandle c);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_value_at(CurveHandle c, double t, [Out] double[] outPoint);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_value_at_many(CurveHandle c, double[] t, int n, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_curvature_at(CurveHandle c, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_curvature_derivative_at(CurveHandle c, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_tangent_at(CurveHandle c, double t, [Out] double[] outVec);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_normal_at(CurveHandle c, double t, [Out] double[] outVec);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_derivative(CurveHandle c, uint n);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_derivative_at(CurveHandle c, uint n, double t, [Out] double[] outVec);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_roots(CurveHandle c, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_extrema(CurveHandle c, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_bounding_box(CurveHandle c, [Out] double[] min, [Out] double[] max);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_split(CurveHandle c, double[] t, int n, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_split_at(CurveHandle c, double t, out IntPtr left, out IntPtr right);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_intersections(CurveHandle c, CurveHandle other, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_project_point(CurveHandle c, double[] point);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_curve_distance(CurveHandle c, double[] point);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_curve_apply_continuity(CurveHandle c, CurveHandle other, double[] beta, int nBeta);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_offset(CurveHandle c, double offset, uint order);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_join(CurveHandle c1, CurveHandle c2, uint order);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_curve_from_polyline(double[] xy, int nPoints, uint order);

    // ============================ PolyCurve ============================
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_new();
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_new_from(IntPtr[] curves, int n);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_copy(PolyCurveHandle p);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_destroy(IntPtr p);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_insert_at(PolyCurveHandle p, uint idx, CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_insert_front(PolyCurveHandle p, CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_insert_back(PolyCurveHandle p, CurveHandle c);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_remove_at(PolyCurveHandle p, uint idx);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_remove_first(PolyCurveHandle p);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_remove_back(PolyCurveHandle p);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern uint bz_polycurve_size(PolyCurveHandle p);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern uint bz_polycurve_curve_idx(PolyCurveHandle p, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_curve(PolyCurveHandle p, uint idx);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_polyline(PolyCurveHandle p, int useFlatness, double flatness, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_polyline_params(PolyCurveHandle p, int useFlatness, double flatness, out int count);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_length(PolyCurveHandle p);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_length_to(PolyCurveHandle p, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_length_between(PolyCurveHandle p, double t1, double t2);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_step(PolyCurveHandle p, double t, double ds);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_end_points(PolyCurveHandle p, [Out] double[] first, [Out] double[] second);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_control_points(PolyCurveHandle p, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_set_control_point(PolyCurveHandle p, uint idx, double[] point);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_value_at(PolyCurveHandle p, double t, [Out] double[] outPoint);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_value_at_many(PolyCurveHandle p, double[] t, int n, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_curvature_at(PolyCurveHandle p, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_curvature_derivative_at(PolyCurveHandle p, double t);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_tangent_at(PolyCurveHandle p, double t, [Out] double[] outVec);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_normal_at(PolyCurveHandle p, double t, [Out] double[] outVec);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_derivative_at(PolyCurveHandle p, uint n, double t, [Out] double[] outVec);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern void bz_polycurve_bounding_box(PolyCurveHandle p, [Out] double[] min, [Out] double[] max);

    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_intersections_curve(PolyCurveHandle p, CurveHandle c, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_intersections_poly(PolyCurveHandle p, PolyCurveHandle other, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_project_point(PolyCurveHandle p, double[] point);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_project_points(PolyCurveHandle p, double[] xy, int nPoints, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern double bz_polycurve_distance(PolyCurveHandle p, double[] point);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_polycurve_distances(PolyCurveHandle p, double[] xy, int nPoints, out int count);

    // ============================ Utils ============================
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_utils_visvalingam_wyatt(double[] xy, int nPoints, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_utils_solve_polynomial(double[] coeffs, int n, out int count);
    [DllImport(Lib, CallingConvention = CallingConvention.Cdecl)]
    internal static extern IntPtr bz_utils_fit_bezier(double[] xy, int nPoints, uint order);
}
