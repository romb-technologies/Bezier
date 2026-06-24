// Copyright 2026 Mirko Kokot. Licensed under the Apache License, Version 2.0.
namespace Bezier;

/// <summary>A point (or vector) in the xy plane. Mirrors Bezier::Point (Eigen::Vector2d).</summary>
public readonly struct Point
{
    public double X { get; }
    public double Y { get; }

    public Point(double x, double y)
    {
        X = x;
        Y = y;
    }

    public override string ToString() => $"({X}, {Y})";
}
