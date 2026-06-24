// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
using Bezier;
using Xunit;

public class CurveTests
{
    // a cubic
    private static Point[] Cp() => new[]
    {
        new Point(0, 0), new Point(1, 2), new Point(3, 3), new Point(4, 0),
    };

    private static void AssertClose(Point a, Point b, double tol = 1e-9)
    {
        Assert.True(System.Math.Abs(a.X - b.X) < tol && System.Math.Abs(a.Y - b.Y) < tol,
            $"expected {b}, got {a}");
    }

    [Fact]
    public void ConstructAndOrder()
    {
        using var c = new Curve(Cp());
        Assert.Equal(3u, c.Order());
        Assert.Equal(4, c.ControlPoints().Length);
    }

    [Fact]
    public void EndpointsMatchValueAt()
    {
        using var c = new Curve(Cp());
        AssertClose(c.ValueAt(0.0), Cp()[0]);
        AssertClose(c.ValueAt(1.0), Cp()[3]);
    }

    [Fact]
    public void BatchValueAtCount()
    {
        using var c = new Curve(Cp());
        Point[] pts = c.ValueAt(new[] { 0.0, 0.5, 1.0 });
        Assert.Equal(3, pts.Length);
    }

    [Fact]
    public void LengthPositive()
    {
        using var c = new Curve(Cp());
        Assert.True(c.Length() > 0.0);
    }

    [Fact]
    public void DerivativeLowersOrder()
    {
        using var c = new Curve(Cp());
        using Curve d1 = c.Derivative();
        using Curve d2 = c.Derivative(2);
        Assert.Equal(2u, d1.Order());
        Assert.Equal(1u, d2.Order());
    }

    [Fact]
    public void SplitRoundtripsEndpoints()
    {
        using var c = new Curve(Cp());
        (Curve left, Curve right) = c.SplitCurve(0.5);
        using (left)
        using (right)
        {
            AssertClose(left.ValueAt(0.0), Cp()[0]);
            AssertClose(right.ValueAt(1.0), Cp()[3]);
            AssertClose(left.ValueAt(1.0), right.ValueAt(0.0)); // C0 at the seam
        }
    }

    [Fact]
    public void BoundingBoxOrdered()
    {
        using var c = new Curve(Cp());
        (Point min, Point max) = c.BoundingBox();
        Assert.True(min.X <= max.X && min.Y <= max.Y);
    }

    [Fact]
    public void PolyCurveSize()
    {
        using var c1 = new Curve(Cp());
        var shifted = System.Array.ConvertAll(Cp(), p => new Point(p.X + 4, p.Y + 4));
        using var c2 = new Curve(shifted);
        using var pc = new PolyCurve(c1, c2);
        Assert.Equal(2u, pc.Size());
        using Curve sub = pc.GetCurve(0);
        Assert.Equal(3u, sub.Order());
    }

    [Fact]
    public void OutOfRangeParamThrowsNotNaN()
    {
        using var c = new Curve(Cp());
        // length(t) with t outside [0,1] throws in C++; the wrapper surfaces it as an exception.
        Assert.Throws<BezierException>(() => c.Length(5.0));
    }
}
